# Test Biostrings conversion functions

test_that("as_DNAString works", {

  dna <- as_DNAString("ACGT")
  expect_s4_class(dna, "DNAString")
  expect_equal(as.character(dna), "ACGT")
})

test_that("as_DNAStringSet works", {

  seqs <- as_DNAStringSet(c("ACGT", "GGCC"))
  expect_s4_class(seqs, "DNAStringSet")
  expect_length(seqs, 2)
})

test_that("as_DNAStringSet accepts names", {

  seqs <- as_DNAStringSet(c("ACGT", "GGCC"), names = c("seq1", "seq2"))
  expect_equal(names(seqs), c("seq1", "seq2"))
})
