#!/bin/sh
set -x

for acc in "$@"; do
  rm -f $acc.png
  apython ../../../pygly/GlycanResource/main.py GlycanImage write $acc
  cp $acc.png ~/www/dropbox/q4dRFkWJuM/image/red-extended-checked
  rm -f $acc.png
  apython ../../../pygly/GlycanResource/main.py GlycanImage write $acc display=compact
  cp $acc.png ~/www/dropbox/q4dRFkWJuM/image/red-compact-checked
  rm -f $acc.png
done
