# BiocRefgetStore (development version)

## New features

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

## Requirements

* Requires `gtars (>= 0.9.2)` -- the release that adds the remote partial-read
  path. With `gtars 0.9.1`, opening a remote store and reading its metadata
  works, but region and whole-sequence reads against a remote store fail.
