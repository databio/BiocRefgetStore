#!/bin/bash
# Rebuild tests/testthat/fixtures/refgenie: the RefgetStore that a refgenie 1.x
# genome folder keeps at <genome_folder>/.refget_store.
#
# Needs refgenie >= 1.0.0a2 on PATH and network access to api.refgenie.org.
# fixtures/rCRSd.fa.gz is refgenie's tests/data/rCRSd.fa, compressed with gzip -n.
# Usage: bash data-raw/refgenie_fixture.sh
set -euo pipefail

PKG_DIR="$(cd "$(dirname "$0")/.." && pwd)"
FIXTURES="$PKG_DIR/tests/testthat/fixtures"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

unset REFGENIE
export REFGENIE_HOME_PATH="$TMP/rghome"
export REFGENIE_SUBSCRIBE_DEFAULT_SERVER=false
export REFGENIE_DISABLE_PLUGINS=all

refgenie init
# Initialized from a FASTA file: the store holds the sequences. Use a copy, since
# reading a FASTA writes an index (.rgsi) beside it.
cp "$FIXTURES/rCRSd.fa.gz" "$TMP/"
refgenie genome init --fasta "$TMP/rCRSd.fa.gz" --name rCRSd --no-build
# Registered from a server without building: the store holds only metadata.
refgenie genome init --server https://api.refgenie.org \
  --digest ia9g4Myony-xxLsSc7eQOC_NdClGgPGu --name R64-1-1 --no-build

rm -rf "$FIXTURES/refgenie"
cp -r "$REFGENIE_HOME_PATH/genomes/.refget_store" "$FIXTURES/refgenie"
refgenie --version
