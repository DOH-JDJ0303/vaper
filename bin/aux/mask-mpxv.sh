#!/bin/bash

INPUT=$1
OUTPUT="$(basename $INPUT | sed 's/.fa.gz//g').masked.fa.gz"
# Mask regions: 1-9400 and 190000-200000
zcat $INPUT | \
awk '
/^>/ {print; next}
{
    seq = $0
    for (i = 1; i <= 9400; i++)  seq = substr(seq, 1, i-1) "N" substr(seq, i+1)
    for (i = 190000; i <= 200000 && i <= length(seq); i++) seq = substr(seq, 1, i-1) "N" substr(seq, i+1)
    print seq
}'  | \
gzip > $OUTPUT

echo "Masked assembly: ${OUTPUT}"