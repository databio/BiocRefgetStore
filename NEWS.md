# BiocRefgetStore (development version)

## New features

* `RefgetGenome.from_refgenie()` opens a genome by name from a refgenie server
  (default `https://api.refgenie.org`). It asks the server which public
  RefgetStore holds the genome and opens it from there. Together with
  `forgeBSgenome()`, this turns any genome the server lists into a BSgenome
  package.

* The docs show how to open a genome from a local refgenie genome folder:
  `RefgetGenome.from_directory("<genome_folder>/.refget_store",
  namespace = "refgenie", alias = ...)`.

## Bug fixes

* `genome[["name"]]` on a sequence whose bases are not in the store (such as a
  refgenie genome registered with `--no-build`) now says so, instead of saying
  the sequence is not in the collection.

## Requirements

* Requires `gtars (>= 0.10.1)`. It refreshes a remote store's manifest each
  time the store is opened; older gtars kept the first cached copy, so genomes
  added to a public store later were not found.

* Imports `jsonlite`.

# BiocRefgetStore 0.1.0

## New features

* `forgeBSgenome()` builds an installable BSgenome data package from a
  `RefgetGenome`, for tools that need a real `BSgenome` object. The collection
  digest is recorded as the package's `genome` field. Requires `BSgenomeForge`.

* `RefgetGenome.from_bsgenome()` loads a BSgenome into a RefgetStore. A genome
  forged with `forgeBSgenome()` and loaded back gets the same collection digest.

* Remote, on-demand sequence retrieval now works end to end against a remote
  `RefgetStore` opened with `RefgetGenome.from_remote()`. Opening fetches only
  collection metadata; each `getSeq()` / `extractRegions()` then pulls just the
  bytes covering the requested region over an HTTP byte-range request -- the
  whole sequence never lands on disk. Region (byte-range) reads are lean and
  persist nothing by default, so each call re-fetches its bytes, which is ideal
  for sparse, random access.

* Coordinate-less whole-sequence access (`genome[["chr"]]` and
  `getSeq(genome, "chr")`) now works against remote stores. When a sequence is
  not resident in a local cache, the full sequence is streamed on demand via
  `gtars::stream_sequence()` instead of returning an empty result.

## Bug fixes

* `getSeq()` on a `RefgetGenome` no longer fails with "unable to find an
  inherited method" after BSgenome or Biostrings is attached. The package now
  adds its method to the Biostrings `getSeq` generic instead of defining its own.

## Requirements

* Biostrings is now a hard dependency (`Imports`). Sequences always come back
  as `DNAString` / `DNAStringSet` unless `as.character = TRUE`.

* Requires `gtars (>= 0.10.0)`.
