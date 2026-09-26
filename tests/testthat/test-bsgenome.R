# Test BSgenome interop: getSeq alongside BSgenome, forging, and the reverse

toy_fasta <- function() {
  fasta_file <- tempfile(fileext = ".fa")
  writeLines(c(
    ">chrA", "ACGTACGTNNACGTAAAACCCCGGGGTTTT",
    ">chrM", "GGGGCCCCAAAATTTTACGT"
  ), fasta_file)
  fasta_file
}

test_that("getSeq still dispatches after BSgenome is attached", {
  skip_if_not_installed("gtars")
  skip_if_not_installed("BSgenome")

  fasta_file <- toy_fasta()
  on.exit(unlink(fasta_file))
  genome <- RefgetGenome.from_fasta(fasta_file)

  suppressPackageStartupMessages(library(BSgenome))
  # Call from the global env, as a user would, so the search-path getSeq is used
  env <- list2env(list(genome = genome), parent = globalenv())
  seq <- eval(quote(getSeq(genome, "chrA", 1, 5)), env)
  expect_equal(as.character(seq), "ACGTA")
})

test_that("forgeBSgenome validates its inputs", {
  skip_if_not_installed("gtars")
  skip_if_not_installed("BSgenomeForge")

  fasta_file <- toy_fasta()
  on.exit(unlink(fasta_file))
  genome <- RefgetGenome.from_fasta(fasta_file)

  expect_error(forgeBSgenome(genome, "Toy", "Toyus example"), "BSgenome\\.")
  expect_error(
    forgeBSgenome(genome, "BSgenome.Toyus.refget.v1", "Toyus example",
                  circ_seqs = "chrZ"),
    "chrZ"
  )
})

test_that("forgeBSgenome and from_bsgenome round-trip to the same digest", {
  skip_if_not_installed("gtars")
  skip_if_not_installed("BSgenomeForge")
  skip_on_cran()

  fasta_file <- toy_fasta()
  dest <- tempfile("forged_")
  lib <- tempfile("forged_lib_")
  dir.create(lib)
  on.exit(unlink(c(fasta_file, dest, lib), recursive = TRUE))
  genome <- RefgetGenome.from_fasta(fasta_file)

  pkg_name <- "BSgenome.Toyus.refget.v1"
  pkg <- forgeBSgenome(genome, pkg_name, organism = "Toyus example",
                       circ_seqs = "chrM", dest = dest)
  expect_true(file.exists(file.path(pkg, "DESCRIPTION")))

  rc <- system2(file.path(R.home("bin"), "R"),
                c("CMD", "INSTALL", "--no-test-load", "-l", shQuote(lib), shQuote(pkg)),
                stdout = FALSE, stderr = FALSE)
  expect_equal(rc, 0L)

  ns <- loadNamespace(pkg_name, lib.loc = lib)
  on.exit(unloadNamespace(pkg_name), add = TRUE)
  bsg <- get("Toyus", envir = ns)

  expect_s4_class(bsg, "BSgenome")
  expect_equal(GenomeInfoDb::seqlengths(bsg), seqlengths(genome))
  expect_equal(unname(GenomeInfoDb::isCircular(bsg)), c(FALSE, TRUE))
  expect_equal(GenomeInfoDb::genome(bsg)[[1]], collection_digest(genome))
  expect_equal(as.character(bsg[["chrA"]]),
               as.character(getSeq(genome, "chrA")))

  back <- RefgetGenome.from_bsgenome(bsg)
  expect_equal(collection_digest(back), collection_digest(genome))
  # The temp FASTA is gone; sequences must still be readable from the store
  expect_equal(getSeq(back, "chrM", as.character = TRUE), "GGGGCCCCAAAATTTTACGT")

  sub <- RefgetGenome.from_bsgenome(bsg, names = "chrA")
  expect_equal(names(sub), "chrA")
})
