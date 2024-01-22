#!/bin/bash

# INPUT="./GSM5067062_2020_6544_human_bedpe.1000lines.bed.gz"
# OUT_DIR="./outputs"

INPUT=$1
OUT_FILE=$2

SAMPLE_NAME="$(basename ${INPUT} .bed.gz)"
echo ${SAMPLE_NAME}

## duplicate rows based on value of column 9
## keep only alignments where column 1 and 4 equal
## swap columns so coordinates are columns 3 to 6
## remove unmapped reads
## get fragment 'start' and 'end' and calculate fragment 'width'
zcat ${INPUT} | \
awk '{ n=$9; while (n--) { print } }' | \
awk -F '\t' '{ if ($1 == $4) { print } }' | \
awk -F '\t' '{ print $1, $4, $2, $3, $5, $6 }' OFS="\t" | \
awk -F '\t' '{ if ($1 != ".") { print } }' | \
awk '{ min=$3;for(i=3;i<=6;i++) if($i<min)min=$i; max=$3;for(i=3;i<=6;i++) if($i>max)max=$i; print $1 "\t0\t0\t" min "\t" max "\t" max-min "\t0" }' \
> ${OUT_FILE}

## EOF
