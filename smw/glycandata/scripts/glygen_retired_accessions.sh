#!/bin/sh

PYTHON=apython
PYGLY=../../../pygly
DATA=../data
EXPORT=../export

set -x

# $PYTHON $PYGLY/GlycanResource/main.py GlyGen allglycans > $DATA/glygen_accessions_alltime.txt
# $PYTHON $PYGLY/GlycanResource/main.py GlyCosmos archived > $DATA/glytoucan_archived.txt
# $PYTHON $PYGLY/GlycanResource/main.py GlyCosmos replaced > $DATA/glytoucan_replaced.txt

$PYTHON ./glygen_retired_accessions.py > $EXPORT/glygen_retired_accessions.txt
