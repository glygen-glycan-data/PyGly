#!/bin/sh
set -x
URL="https://raw.githubusercontent.com/glygener/glycan_dictionary/main/v1.0/glycan_dictionary.csv"
FILE="../data/glycan_dictionary.csv"
rm -f "$FILE"
wget -O "$FILE" "$URL"
