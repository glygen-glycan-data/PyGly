#!/bin/sh

EXPORT="../export"

set -x
./glygen_retired_accessions.py > $EXPORT/glygen_retired_accessions.tsv
