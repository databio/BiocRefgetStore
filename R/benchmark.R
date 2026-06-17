# benchmark.R
#
# Reusable benchmark harness comparing sequence extraction through a
# RefgetStore backend (this package's RefgetGenome) against the standard
# Bioconductor approach. The headline comparator is a REAL BSgenome object
# extracted with BSgenome::getSeq(). By default it is an in-memory BSgenome
# FORGED from the benchmark FASTA (setup_bsgenome_from_fasta), so the bundled
# run always pits RefgetStore against a genuine `is(x, "BSgenome")` object
# with no data-package install and -- because the BSgenome is forged from the
# SAME FASTA the other backends read -- identical coordinates and a meaningful
# cross-backend correctness check. Naming an installed BSgenome data package
# via `bsgenome_pkg` (e.g. "BSgenome.Scerevisiae.UCSC.sacCer3") substitutes
# that real genome instead (point `fasta_path` at its FASTA). An indexed-FASTA
# Rsamtools FaFile is included as an optional secondary baseline.
#
# This file provides reusable functions only -- no plotting, no file I/O.
# A separate vignette consumes the tidy data.frame returned by
# run_benchmarks().
#
### Fairness
#
# The harness is designed so that timing differences reflect the backends'
# extraction machinery, not incidental setup or output-type differences:
#
# * EQUAL COORDINATES. Every backend extracts the SAME regions, produced
#   once by make_benchmark_regions() from seqinfo(genome). Each backend's
#   coordinate convention is handled internally (the RefgetStore path
#   converts 1-based-closed to 0-based-half-open when writing its temp BED;
#   getSeq() on a FaFile/BSgenome with a GRanges is 1-based-closed).
# * EQUAL OUTPUT TYPES. All backends are asked for DNAStringSet/DNAString
#   output (as.character = FALSE). A character-returning call is never
#   timed against a DNAStringSet-returning one.
# * WARM CACHES. Each backend is exercised once (untimed) before timing so
#   first-touch index loads, FFI store loads, and mmap faults are excluded.
#   Reported timings are warm; bench::mark reports the median of N runs.
# * SETUP EXCLUDED. Store construction (RefgetGenome.from_fasta, indexFa,
#   loading the BSgenome) happens in setup_* helpers, outside bench::mark.
# * SAME PROCESS / SAME MACHINE. All contestants run in one R session, so
#   OS, CPU, and disk-cache state are shared.
# * CORRECTNESS CHECKED SEPARATELY. bench::mark(check = FALSE); equality of
#   extracted sequences across backends is asserted ONCE per operation
#   after timing (compared as uppercased character to absorb soft-masking),
#   and recorded in the result's `check` column. A FALSE there signals a
#   coordinate-convention bug and must be investigated.

# -----------------------------------------------------------------------------
# Region construction (shared coordinates for all backends)
# -----------------------------------------------------------------------------

#' Build a shared region set for benchmarking
#'
#' Generates the regions every backend extracts, so all contestants operate
#' on identical 1-based-closed coordinates. Derived from a Seqinfo so it
#' works for both the bundled genome and a real genome.
#'
#' @param seqinfo A `GenomeInfoDb::Seqinfo` (e.g. `seqinfo(genome)`).
#' @param n_single Number of single-region extractions (uses the first).
#' @param n_bulk Number of regions for the bulk extraction test.
#' @param width Region width in bp.
#' @param seed RNG seed for reproducible region placement.
#'
#' @return A list with elements:
#'   `single` (data.frame of `n_single` regions),
#'   `bulk` (data.frame of `n_bulk` regions),
#'   `full_chrom` (character scalar, the longest chromosome name) and
#'   `full_chrom_len` (its length).
#'   Each data.frame has columns `chrom`, `start`, `end` (1-based, closed).
#'
#' @export
make_benchmark_regions <- function(seqinfo, n_single = 1, n_bulk = 1000,
                                   width = 1000, seed = 1) {
  chroms <- GenomeInfoDb::seqnames(seqinfo)
  lens <- GenomeInfoDb::seqlengths(seqinfo)
  names(lens) <- chroms

  # Only use chromosomes long enough to hold a region of `width`.
  usable <- chroms[lens[chroms] > width]
  if (length(usable) == 0) {
    stop("No chromosome is long enough for width = ", width)
  }

  draw <- function(n) {
    ch <- sample(usable, n, replace = TRUE)
    start <- vapply(ch, function(c) {
      sample.int(lens[[c]] - width, 1L)  # 1-based start
    }, integer(1))
    data.frame(chrom = ch, start = start, end = start + width - 1L,
               stringsAsFactors = FALSE)
  }

  withr_seed(seed)
  bulk <- draw(n_bulk)
  single <- draw(n_single)

  full_chrom <- usable[which.max(lens[usable])]

  list(
    single = single,
    bulk = bulk,
    full_chrom = full_chrom,
    full_chrom_len = as.integer(lens[[full_chrom]])
  )
}

# Local seed setter that does not require the withr package.
withr_seed <- function(seed) {
  set.seed(seed)
}

# Build a GRanges from a regions data.frame, attaching seqinfo so it is
# valid against the genome's seqlengths.
.regions_to_granges <- function(df, seqinfo = NULL) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges is required for the benchmark harness")
  }
  # We intentionally do NOT attach seqinfo: getSeq() on a FaFile/BSgenome
  # only needs matching seqnames, and the `seqinfo<-` replacement triggers
  # fragile S4 dispatch when GenomicRanges is loaded via its namespace only
  # (not attached). The `seqinfo` argument is accepted for API symmetry.
  GenomicRanges::GRanges(
    seqnames = df$chrom,
    ranges = IRanges::IRanges(start = df$start, end = df$end)
  )
}

# -----------------------------------------------------------------------------
# Backend setup (all outside the timed loop)
# -----------------------------------------------------------------------------

#' Set up the RefgetStore backend
#'
#' @param fasta_path Path to a FASTA file.
#' @return A `RefgetGenome` object.
#' @export
setup_refget_backend <- function(fasta_path) {
  RefgetGenome.from_fasta(fasta_path)
}

#' Set up the indexed-FASTA (FaFile) baseline backend
#'
#' Ensures a `.fai` index exists (building it once if missing) and returns
#' a `Rsamtools::FaFile`. This is the secondary / fallback baseline; the
#' real BSgenome backend is the primary comparator.
#'
#' @param fasta_path Path to a FASTA file.
#' @return A `Rsamtools::FaFile`, or `NULL` (with a message) if Rsamtools is
#'   unavailable.
#' @export
setup_fafile_backend <- function(fasta_path) {
  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    message("Rsamtools not installed; skipping FaFile baseline.")
    return(NULL)
  }
  fai <- paste0(fasta_path, ".fai")
  if (!file.exists(fai)) {
    # Index a writable copy so we never need to write next to bundled data.
    if (file.access(dirname(fasta_path), 2) != 0) {
      tmp <- tempfile(fileext = ".fa")
      file.copy(fasta_path, tmp, overwrite = TRUE)
      fasta_path <- tmp
    }
    Rsamtools::indexFa(fasta_path)
  }
  Rsamtools::FaFile(fasta_path)
}

#' Set up a real BSgenome backend from an installed data package
#'
#' Loads an installed BSgenome data package by name. The package must be
#' installed separately (it is a large data package and is NOT a Suggests
#' dependency of this package).
#'
#' @param pkgname Name of an installed BSgenome package, e.g.
#'   `"BSgenome.Scerevisiae.UCSC.sacCer3"`.
#' @return The BSgenome object, or `NULL` (with a message) if unavailable.
#' @export
setup_bsgenome_backend <- function(pkgname) {
  if (is.null(pkgname)) return(NULL)
  if (!requireNamespace("BSgenome", quietly = TRUE)) {
    message("BSgenome not installed; skipping BSgenome backend.")
    return(NULL)
  }
  if (!requireNamespace(pkgname, quietly = TRUE)) {
    message(sprintf("BSgenome package '%s' not installed; skipping.", pkgname))
    return(NULL)
  }
  getExportedValue(pkgname, pkgname)
}

#' Forge a real in-memory BSgenome from a FASTA (default primary comparator)
#'
#' Builds a genuine `BSgenome` object from an arbitrary FASTA without
#' installing a data package, by exporting the sequences to a 2bit file in
#' the on-disk layout BSgenome expects (`single_sequences.2bit`) and calling
#' the low-level `BSgenome::BSgenome()` constructor. The returned object is a
#' real `BSgenome` (`is(x, "BSgenome")` is `TRUE`), so `BSgenome::getSeq()`
#' and the `[[` accessor exercise the genuine BSgenome extraction machinery.
#' This is the headline Bioconductor comparator for the bundled run; because
#' it is forged from the SAME FASTA the RefgetStore and FaFile backends read,
#' all backends extract identical coordinates and sequences.
#'
#' @param fasta_path Path to a FASTA file.
#' @return A `BSgenome` object, or `NULL` (with a message) if BSgenome,
#'   rtracklayer, or Biostrings is unavailable, or forging fails.
#' @export
setup_bsgenome_from_fasta <- function(fasta_path) {
  needed <- c("BSgenome", "rtracklayer", "Biostrings")
  have <- vapply(needed, requireNamespace, logical(1), quietly = TRUE)
  if (!all(have)) {
    message(sprintf(
      "Cannot forge BSgenome (%s not installed); skipping forged BSgenome.",
      paste(needed[!have], collapse = ", ")))
    return(NULL)
  }
  tryCatch({
    dna <- Biostrings::readDNAStringSet(fasta_path)
    # FASTA headers may carry descriptions; BSgenome keys on the bare name.
    names(dna) <- sub("\\s.*", "", names(dna))

    seqs_dir <- tempfile("bsgforge")
    dir.create(seqs_dir)
    twobit <- file.path(seqs_dir, "single_sequences.2bit")
    rtracklayer::export.2bit(dna, twobit)

    bsg <- BSgenome::BSgenome(
      organism = "Synthetic benchmark organism",
      common_name = "benchmark",
      provider = "BiocRefgetStore",
      genome = "benchmark",
      release_date = NA,
      source_url = NA,
      seqnames = names(dna),
      circ_seqs = character(0),
      mseqnames = NULL,
      seqs_pkgname = "BSgenome.benchmark",
      seqs_dirpath = seqs_dir
    )
    if (!methods::is(bsg, "BSgenome")) {
      stop("forge did not produce a BSgenome object")
    }
    bsg
  }, error = function(e) {
    message("Forging an in-memory BSgenome failed (", conditionMessage(e),
            "); skipping forged BSgenome backend.")
    NULL
  })
}

# -----------------------------------------------------------------------------
# Per-backend extraction closures (uniform interface, equal output types)
# -----------------------------------------------------------------------------

# Each returns a function of a regions data.frame -> uppercased character
# vector of sequences. Used both for timing and for the correctness check.

.refget_extract_fns <- function(genome, raw = FALSE) {
  if (raw) {
    # Raw extraction: gtars hands back a plain character vector and no
    # Biostrings DNAStringSet is ever constructed. This isolates the actual
    # store-read cost from the downstream Biostrings conversion that the other
    # backends' DNAStringSet output type forces. `as.character = TRUE`
    # short-circuits the DNAStringSet round-trip inside getSeq()/extractRegions().
    return(list(
      single = function(df) {
        getSeq(genome, df$chrom[1], df$start[1], df$end[1], as.character = TRUE)
      },
      bulk = function(df) {
        toupper(extractRegions(genome, df, as.character = TRUE))
      },
      full_chrom = function(chrom) {
        getSeq(genome, chrom, as.character = TRUE)
      }
    ))
  }
  # Converted extraction: build a Biostrings DNAStringSet/DNAString (then back
  # to character), matching the output the BSgenome/FaFile backends return. The
  # delta between this and the raw variant above is exactly the conversion cost.
  list(
    single = function(df) {
      as.character(getSeq(genome, df$chrom[1], df$start[1], df$end[1],
                          as.character = FALSE))
    },
    bulk = function(df) {
      toupper(as.character(extractRegions(genome, df, as.character = FALSE)))
    },
    full_chrom = function(chrom) {
      as.character(getSeq(genome, chrom, as.character = FALSE))
    }
  )
}

.fafile_extract_fns <- function(fafile, seqinfo) {
  list(
    single = function(df) {
      gr <- .regions_to_granges(df[1, , drop = FALSE], seqinfo)
      as.character(Biostrings::getSeq(fafile, gr))
    },
    bulk = function(df) {
      gr <- .regions_to_granges(df, seqinfo)
      toupper(as.character(Biostrings::getSeq(fafile, gr)))
    },
    full_chrom = function(chrom) {
      gr <- GenomicRanges::GRanges(
        seqnames = chrom,
        ranges = IRanges::IRanges(1L, GenomeInfoDb::seqlengths(seqinfo)[[chrom]])
      )
      as.character(Biostrings::getSeq(fafile, gr))
    }
  )
}

.bsgenome_extract_fns <- function(bsg) {
  # getSeq()/[[ on a BSgenome key only on seqnames, so the GRanges needs no
  # attached seqinfo (avoiding a seqinfo,BSgenome dispatch that is fragile
  # when BSgenome is only loaded via its namespace).
  list(
    single = function(df) {
      gr <- .regions_to_granges(df[1, , drop = FALSE])
      toupper(as.character(BSgenome::getSeq(bsg, gr)))
    },
    bulk = function(df) {
      gr <- .regions_to_granges(df)
      toupper(as.character(BSgenome::getSeq(bsg, gr)))
    },
    full_chrom = function(chrom) {
      toupper(as.character(BSgenome::getSeq(bsg, chrom)))
    }
  )
}

# -----------------------------------------------------------------------------
# Core benchmark
# -----------------------------------------------------------------------------

#' Benchmark extraction across backends
#'
#' Times single-region, bulk multi-region, and full-chromosome extraction
#' for each supplied backend over a shared region set, then asserts output
#' correctness across backends (outside the timed loop). Uses `bench::mark`
#' if available, otherwise `microbenchmark`.
#'
#' @param backends A named list of extraction-function lists (as built by
#'   the `setup_*` helpers via [run_benchmarks()]). Each element must have
#'   `single`, `bulk`, and `full_chrom` functions.
#' @param regions Output of [make_benchmark_regions()].
#' @param min_iterations Minimum bench iterations per case (default 10).
#'
#' @return A tidy data.frame with columns `operation`, `backend`,
#'   `n_regions`, `total_bp`, `median_time` (seconds), `mem_alloc` (bytes,
#'   NA without bench), `iterations`, and `check` (logical: output identical
#'   to the reference backend for that operation).
#'
#' @export
benchmark_extraction <- function(backends, regions, min_iterations = 10) {
  use_bench <- requireNamespace("bench", quietly = TRUE)
  use_mb <- !use_bench && requireNamespace("microbenchmark", quietly = TRUE)
  if (!use_bench && !use_mb) {
    stop("Either 'bench' or 'microbenchmark' must be installed.")
  }

  ops <- list(
    single = list(arg = regions$single, n = nrow(regions$single),
                  bp = sum(regions$single$end - regions$single$start + 1L)),
    bulk = list(arg = regions$bulk, n = nrow(regions$bulk),
                bp = sum(regions$bulk$end - regions$bulk$start + 1L)),
    full_chrom = list(arg = regions$full_chrom, n = 1L,
                      bp = regions$full_chrom_len)
  )

  # 1. Warm-up + correctness reference (untimed). Reference is the first
  #    backend; every other backend's output must match it per operation.
  ref_name <- names(backends)[1]
  ref_out <- list()
  checks <- list()  # checks[[op]][[backend]] = logical
  for (op in names(ops)) {
    arg <- ops[[op]]$arg
    for (bk in names(backends)) {
      out <- toupper(backends[[bk]][[op]](arg))  # warm-up call
      if (bk == ref_name) {
        ref_out[[op]] <- out
        checks[[op]][[bk]] <- TRUE
      } else {
        checks[[op]][[bk]] <- isTRUE(all(out == ref_out[[op]]))
      }
    }
  }

  # 2. Timing (check = FALSE so the equality check never contaminates time).
  rows <- list()
  for (op in names(ops)) {
    arg <- ops[[op]]$arg
    for (bk in names(backends)) {
      fn <- backends[[bk]][[op]]
      if (use_bench) {
        res <- bench::mark(
          fn(arg),
          check = FALSE, filter_gc = FALSE,
          min_iterations = min_iterations,
          max_iterations = max(min_iterations, 100)
        )
        median_time <- as.numeric(res$median)
        mem_alloc <- as.numeric(res$mem_alloc)
        iters <- as.integer(res$n_itr)
      } else {
        res <- microbenchmark::microbenchmark(fn(arg), times = min_iterations)
        median_time <- stats::median(res$time) / 1e9  # ns -> s
        mem_alloc <- NA_real_
        iters <- min_iterations
      }
      rows[[length(rows) + 1L]] <- data.frame(
        operation = op,
        backend = bk,
        n_regions = ops[[op]]$n,
        total_bp = ops[[op]]$bp,
        median_time = median_time,
        mem_alloc = mem_alloc,
        iterations = iters,
        check = checks[[op]][[bk]],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# Top-level driver
# -----------------------------------------------------------------------------

#' Run the extraction benchmark
#'
#' Builds the contestants and runs [benchmark_extraction()]. The headline
#' comparator is always a REAL BSgenome object extracted with
#' `BSgenome::getSeq()`:
#'
#' * Default (bundled) run: an in-memory BSgenome is FORGED from the bundled
#'   FASTA (via [setup_bsgenome_from_fasta()]) and used as the primary
#'   comparator. No data-package install and no network are required, and
#'   because the BSgenome is forged from the same FASTA the RefgetStore and
#'   FaFile backends read, all backends extract identical coordinates.
#' * Supplying `bsgenome_pkg`: the named installed BSgenome data package is
#'   used as the primary comparator instead of the forged one (point
#'   `fasta_path` at that genome's FASTA so coordinates match).
#'
#' An indexed-FASTA `FaFile` is added as an optional secondary baseline, and
#' `remote` adds a network-backed RefgetStore backend. If the forged
#' BSgenome cannot be built (e.g. BSgenome/rtracklayer absent) the run does
#' not error: it proceeds with whatever backends are available.
#'
#' @param fasta_path FASTA used to build the RefgetStore, FaFile, and (by
#'   default) forged-BSgenome backends. Defaults to the bundled benchmark
#'   genome.
#' @param bsgenome_pkg Optional name of an installed BSgenome data package to
#'   use as the primary comparator in place of the forged BSgenome. For a
#'   meaningful correctness check, `fasta_path` must be that genome's FASTA.
#' @param remote Optional list with `cache_path`, `remote_url`, and one of
#'   `digest` or (`namespace` + `alias`) to add a `refget_remote` backend.
#' @param n_bulk Number of bulk regions (default 1000).
#' @param width Region width in bp (default 1000).
#' @param min_iterations Minimum bench iterations (default 10).
#'
#' @return The tidy benchmark data.frame from [benchmark_extraction()].
#'
#' @examples
#' \donttest{
#' # Bundled run: RefgetStore vs a real BSgenome forged from the bundled
#' # FASTA (plus a FaFile baseline), no downloads, no data-package install:
#' fa <- system.file("extdata/benchmark/bench_genome.fa",
#'                   package = "BiocRefgetStore")
#' res <- run_benchmarks(fa)
#'
#' # Real-genome run with an installed BSgenome as the primary comparator.
#' # Install the BSgenome package and point fasta_path at its FASTA:
#' # res <- run_benchmarks("sacCer3.fa",
#' #                       bsgenome_pkg = "BSgenome.Scerevisiae.UCSC.sacCer3")
#'
#' # Remote run (adds a network-backed RefgetStore backend):
#' # res <- run_benchmarks(fa, remote = list(
#' #   cache_path = tempfile(),
#' #   remote_url = "https://refgenie.s3.us-east-1.amazonaws.com/pangenome_refget_store",
#' #   digest = "0aHV7I-94paL9Z1H4LNlqsW3WxJhlou5"))
#' }
#'
#' @export
run_benchmarks <- function(fasta_path = system.file("extdata/benchmark/bench_genome.fa",
                                                    package = "BiocRefgetStore"),
                           bsgenome_pkg = NULL, remote = NULL,
                           n_bulk = 1000, width = 1000, min_iterations = 10) {
  if (!nzchar(fasta_path) || !file.exists(fasta_path)) {
    stop("fasta_path does not exist: ", fasta_path)
  }

  # Primary RefgetStore backend + its seqinfo (the shared coordinate space).
  refget <- setup_refget_backend(fasta_path)
  si <- seqinfo(refget)
  regions <- make_benchmark_regions(si, n_bulk = n_bulk, width = width)

  backends <- list()

  # Headline comparator: a REAL BSgenome. Use the named installed data
  # package if supplied, otherwise forge an in-memory BSgenome from the same
  # FASTA so the bundled run always has a genuine BSgenome contestant.
  bsg <- if (!is.null(bsgenome_pkg)) {
    setup_bsgenome_backend(bsgenome_pkg)
  } else {
    setup_bsgenome_from_fasta(fasta_path)
  }
  if (!is.null(bsg)) {
    backends[["bsgenome"]] <- .bsgenome_extract_fns(bsg)
  }

  # Two refget rows: the converted one (DNAStringSet output, apples-to-apples
  # with bsgenome/fafile) and the raw one (gtars character, no conversion), so
  # the table shows how much of refget's time is the Biostrings conversion.
  backends[["refget_converted"]] <- .refget_extract_fns(refget)
  backends[["refget_raw"]]       <- .refget_extract_fns(refget, raw = TRUE)

  # Optional secondary baseline: indexed FASTA via FaFile.
  fafile <- setup_fafile_backend(fasta_path)
  if (!is.null(fafile)) {
    backends[["fafile"]] <- .fafile_extract_fns(fafile, si)
  }

  # Optional remote RefgetStore backend.
  if (!is.null(remote)) {
    rg_remote <- tryCatch(
      RefgetGenome.from_remote(
        cache_path = remote$cache_path,
        remote_url = remote$remote_url,
        digest = remote$digest,
        namespace = remote$namespace,
        alias = remote$alias
      ),
      error = function(e) {
        message("Remote store unavailable (", conditionMessage(e),
                "); skipping refget_remote backend.")
        NULL
      }
    )
    if (!is.null(rg_remote)) {
      backends[["refget_remote"]] <- .refget_extract_fns(rg_remote)
    }
  }

  benchmark_extraction(backends, regions, min_iterations = min_iterations)
}
