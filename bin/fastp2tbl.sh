#!/bin/bash

# fastp2tbl.sh v1.5
# Author: Jared Johnson, jared.johnson@doh.wa.gov

#----- HELP -----#
[ -z $1 ] && echo "fastp2tbl.sh [path/to/raw_fastp.json] [path/to/clean_fastp.json]" && exit

#----- INPUTS -----#
raw_json=$1
clean_json=$2

#----- GATHER METRICS -----#
total_reads_before=$(jq '.summary.before_filtering.total_reads' ${raw_json})
total_bases_before=$(jq '.summary.before_filtering.total_bases' ${raw_json})
read1_mean_length_before=$(jq '.summary.before_filtering.read1_mean_length' ${raw_json})
read2_mean_length_before=$(jq '.summary.before_filtering.read2_mean_length' ${raw_json})
q30_rate_before=$(jq '.summary.before_filtering.q30_rate' ${raw_json})

total_reads_after=$(jq '.summary.before_filtering.total_reads' ${clean_json})
total_bases_after=$(jq '.summary.before_filtering.total_bases' ${clean_json})
read1_mean_length_after=$(jq '.summary.before_filtering.read1_mean_length' ${clean_json})
read2_mean_length_after=$(jq '.summary.before_filtering.read2_mean_length' ${clean_json})
q30_rate_after=$(jq '.summary.before_filtering.q30_rate' ${clean_json})

#----- TABULATE -----#
echo "total_reads_raw,total_bases_raw,read1_mean_length_raw,read2_mean_length_raw,q30_rate_raw,total_reads_clean,total_bases_clean,read1_mean_length_clean,read2_mean_length_clean,q30_rate_clean"
echo "${total_reads_before},${total_bases_before},${read1_mean_length_before},${read2_mean_length_before},${q30_rate_before},${total_reads_after},${total_bases_after},${read1_mean_length_after},${read2_mean_length_after},${q30_rate_after}"
