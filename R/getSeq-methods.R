# Define getSeq generic if not already available from Biostrings/BSgenome
# This allows the package to work without Biostrings installed
if (!isGeneric("getSeq")) {
  setGeneric("getSeq", function(x, ...) standardGeneric("getSeq"))
}

#' getSeq method for RefgetGenome
#'
#' Extract sequences from a RefgetGenome object using BSgenome-compatible syntax.
#'
#' @param x A RefgetGenome object
#' @param names Sequence names (character vector) or a GRanges object
#' @param start Start positions (integer vector, 1-based)
#' @param end End positions (integer vector, 1-based)
#' @param strand Strand ("+" or "-"). Default is "+".
#' @param as.character If TRUE, return character strings instead of DNAString/DNAStringSet
#' @param ... Additional arguments (ignored)
#'
#' @return
#' - Single sequence: DNAString (or character if as.character=TRUE or Biostrings unavailable)
#' - Multiple sequences: DNAStringSet (or character vector)
#'
#' @examples
#' \dontrun{
#' genome <- RefgetGenome.from_fasta("genome.fa")
#'
#' # Full chromosome
#' seq <- getSeq(genome, "chr1")
#'
#' # Region by coordinates
#' seq <- getSeq(genome, "chr1", 1000, 2000)
#'
#' # Multiple regions
#' seqs <- getSeq(genome, c("chr1", "chr2"), c(1000, 5000), c(2000, 6000))
#'
#' # From GRanges (requires GenomicRanges)
#' library(GenomicRanges)
#' gr <- GRanges(c("chr1:1000-2000", "chr2:5000-6000:-"))
#' seqs <- getSeq(genome, gr)
#' }
#'
#' @rdname getSeq-methods
#' @aliases getSeq
#' @export
setMethod("getSeq", "RefgetGenome",
  function(x, names, start = NA, end = NA, strand = "+", as.character = FALSE, ...) {
    # Handle GRanges input
    if (inherits(names, "GRanges")) {
      return(.getSeq_GRanges(x, names, as.character = as.character))
    }

    # Ensure names is character
    names <- as.character(names)

    # Handle single sequence (full or region)
    if (length(names) == 1 && length(start) <= 1 && length(end) <= 1) {
      return(.getSeq_single(x, names, start, end, strand, as.character))
    }

    # Vectorized extraction
    .getSeq_vectorized(x, names, start, end, strand, as.character)
  }
)

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Extract a single sequence or region
#' @keywords internal
.getSeq_single <- function(genome, name, start, end, strand, as.character) {
  if (!is.na(start) && !is.na(end)) {
    # Coordinate extraction: indexed substring via gtars, NOT a whole-sequence
    # fetch. gtars uses 0-based half-open (BED) coordinates, so start - 1.
    retrieved <- gtars::get_seqs_from_regions(
      genome@store, genome@collection_digest,
      name, as.integer(start - 1L), as.integer(end)
    )
    seq_string <- retrieved$sequence[[1]]
  } else {
    # Whole-sequence request (no coordinates): fetch the full sequence.
    record <- gtars::get_sequence_by_name(genome@store, genome@collection_digest, name)
    if (is.null(record)) {
      stop(sprintf("Sequence '%s' not found in collection", name))
    }
    seq_string <- record@data
  }

  # Handle negative strand
  if (identical(strand, "-")) {
    seq_string <- .reverse_complement(seq_string)
  }

  # Convert to DNAString if requested and Biostrings available
  if (!as.character && requireNamespace("Biostrings", quietly = TRUE)) {
    return(Biostrings::DNAString(seq_string))
  }

  seq_string
}

#' Vectorized sequence extraction
#' @keywords internal
.getSeq_vectorized <- function(genome, names, start, end, strand, as.character) {
  n <- length(names)

  # Recycle start/end/strand to match names length
  if (length(start) == 1) start <- rep(start, n)
  if (length(end) == 1) end <- rep(end, n)
  if (length(strand) == 1) strand <- rep(strand, n)

  if (length(start) != n || length(end) != n || length(strand) != n) {
    stop("Length mismatch: names, start, end, and strand must have compatible lengths")
  }

  # Use in-memory bulk extraction if all regions have coordinates
  if (!any(is.na(start)) && !any(is.na(end))) {
    # Indexed extraction via a single FFI call: region vectors in, flat
    # character vector out (no temp BED file, no per-region object boxing).
    # gtars uses 0-based half-open (BED) coordinates, so start - 1.
    retrieved <- gtars::get_seqs_from_regions(
      genome@store, genome@collection_digest,
      names, as.integer(start - 1L), as.integer(end)
    )

    seqs <- retrieved$sequence

    # Handle negative strand
    minus_idx <- which(strand == "-")
    if (length(minus_idx) > 0) {
      seqs[minus_idx] <- vapply(seqs[minus_idx], .reverse_complement, character(1))
    }

    result_names <- sprintf("%s:%d-%d", names, start, end)
  } else {
    # Fallback: full chromosome extraction (no coordinates)
    seqs <- vapply(seq_len(n), function(i) {
      .getSeq_single(genome, names[i], start[i], end[i], strand[i], as.character = TRUE)
    }, character(1))

    result_names <- names
  }

  names(seqs) <- result_names

  # Convert to DNAStringSet if requested
  if (!as.character && requireNamespace("Biostrings", quietly = TRUE)) {
    return(Biostrings::DNAStringSet(seqs))
  }

  seqs
}

#' Extract sequences from GRanges
#' @keywords internal
.getSeq_GRanges <- function(genome, gr, as.character) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges package required for GRanges input")
  }

  # Extract components from GRanges
  names <- as.character(GenomicRanges::seqnames(gr))
  start <- GenomicRanges::start(gr)
  end <- GenomicRanges::end(gr)
  strand <- as.character(GenomicRanges::strand(gr))
  strand[strand == "*"] <- "+"  # Treat unstranded as +

  .getSeq_vectorized(genome, names, start, end, strand, as.character)
}
