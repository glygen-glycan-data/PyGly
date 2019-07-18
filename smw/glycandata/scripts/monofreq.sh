#!/bin/sh

export LC_ALL=C
awk 'NF == 3 {print $2,$3}' | sort | uniq -c | awk '{print $3,$2,$1}' | sort -k1,1 -k3rn,3 -k2,2 | tr ' ' '\t' 
