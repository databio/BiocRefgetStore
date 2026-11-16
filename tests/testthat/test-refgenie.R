# Genomes managed by refgenie (>= 1.0): a local genome folder, and a refgenie
# server. fixtures/refgenie is a store written by refgenie itself; rebuild it
# with data-raw/refgenie_fixture.sh.

RCRSD_DIGEST <- "jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP"
YEAST_DIGEST <- "ia9g4Myony-xxLsSc7eQOC_NdClGgPGu"

# Copy the fixture to <tmp>/genomes/.refget_store, the layout of a genome
# folder, so tests never touch the installed package.
refgenie_store <- function() {
  genomes <- file.path(tempfile("refgenie_"), "genomes")
  dir.create(genomes, recursive = TRUE)
  file.copy(test_path("fixtures", "refgenie"), genomes, recursive = TRUE)
  file.rename(file.path(genomes, "refgenie"), file.path(genomes, ".refget_store"))
  file.path(genomes, ".refget_store")
}

rcrsd_fasta <- function() {
  fasta <- tempfile(fileext = ".fa.gz")
  file.copy(test_path("fixtures", "rCRSd.fa.gz"), fasta)
  fasta
}

# Install a forged package into a temp library and return its BSgenome object.
install_forged <- function(pkg, pkg_name, lib) {
  rc <- system2(file.path(R.home("bin"), "R"),
                c("CMD", "INSTALL", "--no-test-load", "-l", shQuote(lib), shQuote(pkg)),
                stdout = FALSE, stderr = FALSE)
  expect_equal(rc, 0L)
  ns <- loadNamespace(pkg_name, lib.loc = lib)
  get(strsplit(pkg_name, ".", fixed = TRUE)[[1]][2], envir = ns)
}

skip_if_no_refgenie_server <- function() {
  skip_if_not_installed("gtars")
  skip_on_cran()
  if (!nzchar(Sys.getenv("BIOCREFGETSTORE_RUN_REMOTE"))) {
    skip("Set BIOCREFGETSTORE_RUN_REMOTE=1 to run live refgenie server tests")
  }
  skip_if_offline("api.refgenie.org")
}

test_that("a refgenie genome folder opens by refgenie name", {
  skip_if_not_installed("gtars")
  path <- refgenie_store()
  on.exit(unlink(dirname(dirname(path)), recursive = TRUE))

  genome <- RefgetGenome.from_directory(path, namespace = "refgenie", alias = "rCRSd")
  expect_equal(collection_digest(genome), RCRSD_DIGEST)
  expect_equal(seqlengths(genome), c(rCRSd = 33138L))

  fasta <- rcrsd_fasta()
  on.exit(unlink(fasta), add = TRUE)
  ref <- Biostrings::readDNAStringSet(fasta)[[1]]
  expect_equal(as.character(genome[["rCRSd"]]), as.character(ref))
  expect_equal(as.character(getSeq(genome, "rCRSd", 101, 160)),
               as.character(Biostrings::subseq(ref, 101, 160)))

  # The refgenie genome digest is the seqcol digest of the FASTA
  expect_equal(collection_digest(RefgetGenome.from_fasta(fasta)), RCRSD_DIGEST)
})

test_that("refgenie names are not in the 'name' namespace of a genome folder", {
  skip_if_not_installed("gtars")
  path <- refgenie_store()
  on.exit(unlink(dirname(dirname(path)), recursive = TRUE))

  expect_error(RefgetGenome.from_directory(path, namespace = "name", alias = "rCRSd"),
               "No collection found")
})

test_that("a genome registered without its sequences gives a clear error", {
  skip_if_not_installed("gtars")
  path <- refgenie_store()
  on.exit(unlink(dirname(dirname(path)), recursive = TRUE))

  genome <- RefgetGenome.from_directory(path, namespace = "refgenie", alias = "R64-1-1")
  expect_equal(collection_digest(genome), YEAST_DIGEST)
  expect_equal(length(genome), 17L)
  expect_error(suppressMessages(genome[["Mito"]]), "could not be loaded from this store")
  expect_error(genome[["chrNope"]], "not found in collection")
})

test_that("forgeBSgenome works on a genome from a refgenie folder", {
  skip_if_not_installed("gtars")
  skip_if_not_installed("BSgenomeForge")
  skip_on_cran()
  path <- refgenie_store()
  dest <- tempfile("forged_")
  lib <- tempfile("forged_lib_")
  dir.create(lib)
  on.exit(unlink(c(dirname(dirname(path)), dest, lib), recursive = TRUE))

  genome <- RefgetGenome.from_directory(path, namespace = "refgenie", alias = "rCRSd")
  pkg_name <- "BSgenome.Hsapiens.refgenie.rCRSd"
  pkg <- forgeBSgenome(genome, pkg_name, organism = "Homo sapiens",
                       circ_seqs = "rCRSd", dest = dest)
  bsg <- install_forged(pkg, pkg_name, lib)
  on.exit(unloadNamespace(pkg_name), add = TRUE)

  expect_equal(GenomeInfoDb::genome(bsg)[[1]], RCRSD_DIGEST)
  expect_equal(as.character(bsg[["rCRSd"]]), as.character(genome[["rCRSd"]]))
})

test_that("RefgetGenome.from_refgenie opens a genome by name from the server", {
  skip_if_no_refgenie_server()
  cache <- tempfile("refgenie_cache_")
  on.exit(unlink(cache, recursive = TRUE))

  genome <- RefgetGenome.from_refgenie("rCRSd", cache_path = cache)
  expect_equal(collection_digest(genome), RCRSD_DIGEST)

  fasta <- rcrsd_fasta()
  on.exit(unlink(fasta), add = TRUE)
  ref <- Biostrings::readDNAStringSet(fasta)[[1]]
  expect_equal(as.character(genome[["rCRSd"]]), as.character(ref))

  by_digest <- RefgetGenome.from_refgenie(RCRSD_DIGEST, cache_path = cache)
  expect_equal(collection_digest(by_digest), RCRSD_DIGEST)

  # A name that is 32 characters long, like a digest, is still looked up as a name
  long_name <- "vigna_radiata_var_radiata_LIS_v6"
  by_long_name <- RefgetGenome.from_refgenie(long_name, cache_path = file.path(cache, "other"))
  expect_false(collection_digest(by_long_name) == long_name)

  expect_error(RefgetGenome.from_refgenie("no-such-genome-xyz", cache_path = cache),
               "no-such-genome-xyz")
})

test_that("a BSgenome package forges straight from the refgenie server", {
  skip_if_no_refgenie_server()
  skip_if_not_installed("BSgenomeForge")
  cache <- tempfile("refgenie_cache_")
  dest <- tempfile("forged_")
  lib <- tempfile("forged_lib_")
  dir.create(lib)
  on.exit(unlink(c(cache, dest, lib), recursive = TRUE))

  genome <- RefgetGenome.from_refgenie("rCRSd", cache_path = cache)
  pkg_name <- "BSgenome.Hsapiens.refgenie.rCRSdremote"
  pkg <- forgeBSgenome(genome, pkg_name, organism = "Homo sapiens",
                       circ_seqs = "rCRSd", dest = dest)
  bsg <- install_forged(pkg, pkg_name, lib)
  on.exit(unloadNamespace(pkg_name), add = TRUE)

  expect_equal(GenomeInfoDb::genome(bsg)[[1]], RCRSD_DIGEST)
  expect_equal(collection_digest(RefgetGenome.from_bsgenome(bsg)), RCRSD_DIGEST)
})
