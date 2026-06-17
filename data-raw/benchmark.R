#!/usr/bin/env Rscript
# Regenerate the committed benchmark artifacts on the full human genome (GRCh38):
#   vignettes/benchmark_results.rds   (tidy results data frame)
#   vignettes/benchmark_plot.png      (faceted comparison plot)
#
# Run from the package root:
#   Rscript data-raw/benchmark.R
#
# Requires the benchmark Suggests (bench, BSgenome, rtracklayer, Biostrings,
# ggplot2) and the BSgenome.Hsapiens.UCSC.hg38 data package. Exports the FASTA to
# tempdir() (a few GB of scratch, removed at the end).
#
# To benchmark a different assembly, change `bsgenome_name` / `chroms`.

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)   # use the working-tree code
  } else {
    library(BiocRefgetStore)
  }
  library(Biostrings)
  library(ggplot2)
})

bsgenome_name <- "BSgenome.Hsapiens.UCSC.hg38"
chroms        <- paste0("chr", c(1:22, "X", "Y", "M"))   # primary assembly
n_bulk        <- 2000

suppressMessages(library(bsgenome_name, character.only = TRUE))
bsg <- get(bsgenome_name)

scratch <- file.path(tempdir(), "biocrefgetstore_bench")
dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
fa <- file.path(scratch, "genome.fa")

# Export the BSgenome's own sequence to FASTA one chromosome at a time (append),
# so the RefgetStore and the BSgenome comparator are byte-identical (the
# correctness check is meaningful) and peak memory stays ~one chromosome.
if (!file.exists(fa)) {
  message("Exporting ", length(chroms), " chromosomes -> ", fa)
  for (i in seq_along(chroms)) {
    cn <- chroms[i]
    ss <- DNAStringSet(bsg[[cn]]); names(ss) <- cn
    writeXStringSet(ss, fa, append = (i > 1))
    rm(ss); gc(verbose = FALSE)
  }
}

res <- run_benchmarks(
  fasta_path     = fa,
  bsgenome_pkg   = bsgenome_name,
  n_bulk         = n_bulk,
  width          = 1000,
  min_iterations = 5
)

saveRDS(res, "vignettes/benchmark_results.rds")

p <- ggplot(res, aes(x = backend, y = median_time * 1e3, fill = backend)) +
  geom_col() +
  facet_wrap(~ operation, scales = "free_y") +
  labs(x = NULL, y = "Median time (ms)",
       title = "Extraction time by backend, faceted by operation (GRCh38/hg38)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("vignettes/benchmark_plot.png", p, width = 7, height = 4, dpi = 120)

print(res[, c("operation", "backend", "n_regions", "total_bp", "median_time", "check")])
unlink(scratch, recursive = TRUE)
message("Wrote vignettes/benchmark_results.rds and vignettes/benchmark_plot.png")
