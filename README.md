# BiocRefgetStore

BSgenome-compatible interface to RefgetStore for R/Bioconductor.

## Overview

BiocRefgetStore provides a bridge between the gtars RefgetStore format and Bioconductor's BSgenome API. This allows you to use RefgetStore-backed genomes with the familiar `getSeq()` interface that Bioconductor users expect.

## Installation

```r
# Install gtars first (required)
# See: https://github.com/databio/gtars

# Install BiocRefgetStore
devtools::install_github("databio/BiocRefgetStore")
# Or from local source:
# devtools::install_local("path/to/BiocRefgetStore")
```

## Quick Start

```r
library(BiocRefgetStore)

# Create a genome from a FASTA file
genome <- RefgetGenome.from_fasta("genome.fa")

# BSgenome-compatible access
seq <- getSeq(genome, "chr1", 1000, 2000)           # Returns DNAString
seqs <- getSeq(genome, c("chr1", "chr2"))           # Returns DNAStringSet
genome[["chr1"]]                                    # Full chromosome

# Standard accessors
seqnames(genome)         # c("chr1", "chr2", ...)
seqlengths(genome)       # Named integer vector
seqinfo(genome)          # Seqinfo object
```

## Loading Genomes

```r
# From FASTA file (creates in-memory store)
genome <- RefgetGenome.from_fasta("/path/to/genome.fa")

# From persisted RefgetStore directory
genome <- RefgetGenome.from_directory("/path/to/store", digest = "abc123...")

# From RefgetStore with alias
store <- gtars::refget_store_open_local("/path/to/store")
genome <- RefgetGenome(store, namespace = "refseq", alias = "GRCh38")

# From a remote store (cloud-backed, with on-demand local caching)
genome <- RefgetGenome.from_remote(
  cache_path = "~/.cache/refget/pangenome",
  remote_url = "https://refgenie.s3.us-east-1.amazonaws.com/pangenome_refget_store",
  digest = "0qveCdMlbF_kYn6XWb7YBy-FtRZ6gSAL"
)
```

## Remote store, no local files

You can work entirely against a remote RefgetStore: point at its URL plus a
local cache directory and pull only the bytes you touch. Opening the store
fetches just metadata; each `getSeq()` downloads only the sequence data it
needs, and that data is cached under `cache_path` so later calls (and later
sessions with the same `cache_path`) reuse it instead of re-downloading.

```r
library(BiocRefgetStore)

pangenome_url <- "https://refgenie.s3.us-east-1.amazonaws.com/pangenome_refget_store"

# 1. Open: only metadata is fetched here -- no sequence data, no local store.
genome <- RefgetGenome.from_remote(
  cache_path = "~/.cache/refget/pangenome",
  remote_url = pangenome_url,
  digest = "0qveCdMlbF_kYn6XWb7YBy-FtRZ6gSAL"
)
genome
seqnames(genome)[1:5]
seqlengths(genome)[1:5]

# Pick a real sequence name from the collection.
chrom <- seqnames(genome)[1]

# 2. Single region -- pulls just this slice, then caches it.
getSeq(genome, chrom, start = 1000, end = 2000)

# 3. Vectorized extraction across several regions in one call.
getSeq(genome,
       names = rep(chrom, 3),
       start = c(100, 5000, 20000),
       end   = c(199, 5099, 20099))

# 4. Full sequence (use a small one to bound the download).
small <- names(which.min(seqlengths(genome)))
genome[[small]]

# 5. Bulk extraction from a data.frame of regions.
regions <- data.frame(
  chrom = rep(chrom, 2),
  start = c(100, 5000),
  end   = c(199, 5099)
)
extractRegions(genome, regions)
extractToFasta(genome, regions, "regions.fa")
```

> **Note (gtars 0.9.0):** opening remote stores and reading their metadata
> works today, but on-demand sequence-*data* retrieval against this public
> store is currently blocked by a gtars 0.9.0 bug (it does not strip the
> `SQ.` digest prefix when building the remote `.seq` object path, so data
> fetches 403). The data steps above are the intended workflow and become
> functional once that prefix handling is fixed upstream in gtars.

For the same workflow in Python, see
[`data_loaders/demo_remote_store.py`](https://github.com/refgenie/refget/blob/dev/data_loaders/demo_remote_store.py)
in the `refgenie/refget` repository.

## Sequence Extraction

```r
# Single region
seq <- getSeq(genome, "chr1", 1000, 2000)

# Multiple regions
seqs <- getSeq(genome,
               names = c("chr1", "chr2"),
               start = c(1000, 5000),
               end = c(2000, 6000))

# From GRanges (requires GenomicRanges)
library(GenomicRanges)
gr <- GRanges(c("chr1:1000-2000", "chr2:5000-6000:-"))
seqs <- getSeq(genome, gr)

# Strand-aware extraction
seq <- getSeq(genome, "chr1", 1000, 2000, strand = "-")  # Reverse complement
```

## Bulk Extraction

```r
# Extract multiple regions efficiently
regions <- data.frame(
  chrom = c("chr1", "chr2", "chr3"),
  start = c(1000, 5000, 10000),
  end = c(2000, 6000, 11000)
)
seqs <- extractRegions(genome, regions)

# Write regions to FASTA
extractToFasta(genome, regions, "extracted.fa")

# Export specific chromosomes
exportChromosomes(genome, c("chr1", "chr2"), "subset.fa")
```

## RefgetStore-Specific Features

```r
# Get the seqcol digest
collection_digest(genome)

# Get coordinate system (for compatibility checking)
coordinate_system(genome)

# Get per-sequence digests
sequence_digests(genome)

# Access underlying RefgetStore
store <- store(genome)
```

## Dependencies

- **Required**: gtars (for RefgetStore), GenomeInfoDb (for Seqinfo)
- **Optional**: Biostrings (for DNAString/DNAStringSet), GenomicRanges (for GRanges support)

Without Biostrings, sequences are returned as character strings.
