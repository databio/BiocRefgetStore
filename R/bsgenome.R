#' Forge a BSgenome data package from a RefgetGenome
#'
#' Builds an installable BSgenome data package (the kind named like
#' \code{BSgenome.Hsapiens.UCSC.hg38}) from the sequences in a RefgetGenome.
#' Use this when a tool needs a real \code{BSgenome} object; for plain sequence
#' access, a RefgetGenome already supports \code{getSeq()} directly.
#'
#' Each sequence is exported from the RefgetStore to a temporary FASTA file and
#' handed to \code{BSgenomeForge::forgeBSgenomeDataPkg()}, which packs them into
#' a single \code{.2bit} file inside the new package. The collection digest is
#' recorded as the package's \code{genome} field, so the forged package can be
#' traced back to the exact RefgetStore collection it came from.
#'
#' @param genome A RefgetGenome object
#' @param pkg_name Package name, following the BSgenome convention
#'   \code{BSgenome.<Abbrev>.<provider>.<build>} (e.g.
#'   \code{"BSgenome.Hsapiens.refget.GRCh38"}). The second field becomes the
#'   name of the BSgenome object the package exports.
#' @param organism Scientific name of the organism (e.g. \code{"Homo sapiens"})
#' @param common_name Common name of the organism (default: \code{organism})
#' @param circ_seqs Names of circular sequences such as \code{"chrM"}
#'   (default: none)
#' @param version Version of the forged package (default: \code{"1.0.0"})
#' @param dest Directory in which to create the package (default: \code{"."})
#' @param replace If TRUE, overwrite an existing package directory
#'
#' @return Invisibly, the path to the new package directory. Install it with
#'   \code{R CMD INSTALL <path>}.
#'
#' @seealso \code{\link{RefgetGenome.from_bsgenome}} for the reverse direction
#'
#' @examples
#' \dontrun{
#' genome <- RefgetGenome.from_fasta("genome.fa")
#' pkg <- forgeBSgenome(genome, "BSgenome.Hsapiens.refget.GRCh38",
#'                      organism = "Homo sapiens", circ_seqs = "chrM",
#'                      dest = tempdir())
#' install.packages(pkg, repos = NULL, type = "source")
#' }
#'
#' @export
forgeBSgenome <- function(genome, pkg_name, organism, common_name = organism,
                          circ_seqs = character(0), version = "1.0.0",
                          dest = ".", replace = FALSE) {
  if (!requireNamespace("BSgenomeForge", quietly = TRUE)) {
    stop("The 'BSgenomeForge' package is required for forgeBSgenome(). ",
         "Install with: BiocManager::install('BSgenomeForge')", call. = FALSE)
  }
  pkg_fields <- strsplit(pkg_name, ".", fixed = TRUE)[[1]]
  if (length(pkg_fields) != 4 || pkg_fields[1] != "BSgenome") {
    stop("pkg_name must look like 'BSgenome.<Abbrev>.<provider>.<build>'")
  }
  seq_names <- names(genome)
  bad_circ <- setdiff(circ_seqs, seq_names)
  if (length(bad_circ) > 0) {
    stop("circ_seqs not found in genome: ", paste(bad_circ, collapse = ", "))
  }
  digest <- collection_digest(genome)

  # One FASTA per sequence, named <seqname>.fa, as BSgenomeForge expects
  src_dir <- tempfile("forge_src_")
  dir.create(src_dir)
  on.exit(unlink(src_dir, recursive = TRUE), add = TRUE)
  for (nm in seq_names) {
    exportChromosomes(genome, names = nm,
                      output_path = file.path(src_dir, paste0(nm, ".fa")))
  }

  # seqnames and circ_seqs are required: the build is a refget digest, which
  # GenomeInfoDb doesn't know, so BSgenomeForge can't look them up itself.
  seed <- list(
    Package = pkg_name,
    Title = sprintf("Full genome sequences for %s from a RefgetStore", organism),
    Description = sprintf(
      "Full genome sequences for %s, forged from RefgetStore collection %s.",
      organism, digest),
    Version = version,
    organism = organism,
    common_name = common_name,
    provider = "refget",
    genome = digest,
    release_date = format(Sys.Date(), "%Y-%m-%d"),
    source_url = "https://ga4gh.github.io/refget/",
    organism_biocview = gsub(" ", "_", organism),
    BSgenomeObjname = pkg_fields[2],
    seqnames = deparse1(seq_names),
    circ_seqs = deparse1(as.character(circ_seqs)),
    seqs_srcdir = src_dir
  )
  seed_file <- file.path(src_dir, "seed.dcf")
  # No line wrapping: BSgenomeForge pastes values into its DESCRIPTION template
  # verbatim, so a wrapped continuation line would corrupt it.
  write.dcf(as.data.frame(seed, check.names = FALSE), seed_file, width = 10000)

  if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)
  BSgenomeForge::forgeBSgenomeDataPkg(seed_file, destdir = dest,
                                      replace = replace, verbose = FALSE)
  invisible(file.path(dest, pkg_name))
}

#' Create RefgetGenome from a BSgenome
#'
#' Loads the sequences of a BSgenome object into a new in-memory RefgetStore,
#' giving it a refget collection digest. The reverse of
#' \code{\link{forgeBSgenome}}: a genome forged from a RefgetGenome and loaded
#' back this way gets the same collection digest.
#'
#' Sequences are written one at a time to a temporary FASTA file, so only one
#' sequence is held in R memory at once.
#'
#' @param bsgenome A BSgenome object (e.g. \code{BSgenome.Hsapiens.UCSC.hg38})
#' @param names Optional sequence names to include (NULL = all)
#'
#' @return A RefgetGenome object
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' genome <- RefgetGenome.from_bsgenome(BSgenome.Hsapiens.UCSC.hg38)
#' collection_digest(genome)
#' }
#'
#' @export
RefgetGenome.from_bsgenome <- function(bsgenome, names = NULL) {
  if (!requireNamespace("BSgenome", quietly = TRUE)) {
    stop("The 'BSgenome' package is required for RefgetGenome.from_bsgenome(). ",
         "Install with: BiocManager::install('BSgenome')", call. = FALSE)
  }
  if (!is(bsgenome, "BSgenome")) {
    stop("bsgenome must be a BSgenome object")
  }
  names <- names %||% GenomeInfoDb::seqnames(bsgenome)

  fasta <- tempfile(fileext = ".fa")
  on.exit(unlink(fasta), add = TRUE)
  for (i in seq_along(names)) {
    seq <- Biostrings::DNAStringSet(getSeq(bsgenome, names[i]))
    base::names(seq) <- names[i]
    Biostrings::writeXStringSet(seq, fasta, append = i > 1)
  }
  RefgetGenome.from_fasta(fasta)
}
