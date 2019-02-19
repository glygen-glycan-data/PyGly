#!/bin/sh

awk '{print $2}' | sort -u | $HOME/projects/ibpep/peptide_scan -i $2 -P - 
