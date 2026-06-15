#!/usr/bin/env Rscript
# make_benchmark_fasta.R
#
# Deterministically generate the small benchmark genome bundled in
# inst/extdata/benchmark/. The generated FASTA is committed so the
# benchmark/vignette never needs to regenerate it at build time; this
# script exists only so the data is reproducible/regenerable.
#
# Regenerate with:
#   Rscript inst/scripts/make_benchmark_fasta.R
# (run from the package root). Writes:
#   inst/extdata/benchmark/bench_genome.fa
#   inst/extdata/benchmark/bench_genome.fa.fai   (if Rsamtools available)
#   inst/extdata/benchmark/bench_regions.bed
#
# Design: a few chromosomes totaling a few Mb -- large enough that
# full-chromosome and bulk extraction take measurable time, small enough
# to commit and keep the build fast. All output is deterministic (fixed
# RNG seed) and uppercase ACGT (no soft-masking), so correctness checks
# across backends compare cleanly.

set.seed(20260615L)

# Chromosome lengths (bp). Total ~3.5 Mb.
chrom_lengths <- c(chr1 = 2000000L, chr2 = 1000000L, chr3 = 500000L)

bases <- c("A", "C", "G", "T")

random_seq <- function(n) {
  paste(sample(bases, n, replace = TRUE), collapse = "")
}

# --- Write FASTA (80 bp lines) ---
out_dir <- file.path("inst", "extdata", "benchmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fa_path <- file.path(out_dir, "bench_genome.fa")

con <- file(fa_path, open = "wt")
line_width <- 80L
for (chrom in names(chrom_lengths)) {
  len <- chrom_lengths[[chrom]]
  seq <- random_seq(len)
  writeLines(paste0(">", chrom), con)
  starts <- seq.int(1L, len, by = line_width)
  lines <- substring(seq, starts, pmin(starts + line_width - 1L, len))
  writeLines(lines, con)
}
close(con)
cat("Wrote", fa_path, "\n")

# --- Write FAI index (so the FaFile baseline needs no write step) ---
if (requireNamespace("Rsamtools", quietly = TRUE)) {
  Rsamtools::indexFa(fa_path)
  cat("Wrote", paste0(fa_path, ".fai"), "\n")
} else {
  cat("Rsamtools not available; skipping .fai (harness builds it on demand)\n")
}

# --- Write a deterministic BED of bulk benchmark regions ---
# 1000 regions of width 1000, spread across all chromosomes.
n_bulk <- 1000L
width <- 1000L
bed_path <- file.path(out_dir, "bench_regions.bed")

chroms <- sample(names(chrom_lengths), n_bulk, replace = TRUE)
starts0 <- vapply(chroms, function(ch) {
  sample.int(chrom_lengths[[ch]] - width, 1L) - 1L  # 0-based
}, integer(1))
ends <- starts0 + width

bed <- data.frame(chrom = chroms, start = starts0, end = ends)
write.table(bed, bed_path, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
cat("Wrote", bed_path, "(", n_bulk, "regions )\n")
