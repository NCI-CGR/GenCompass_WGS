#!/bin/bash

# Usage: ./prepend_INFO_labels.sh <vcf_without_headers> <vcf_header_file> <caller_label>
#
# Takes a VCF, uncompressed and split into vars and header files,
# and returns, uncompressed:
#   1. a variant file where all input fields (;-delimited)
#   and format fields (:-delimited) are prepended with the caller
#   label (e.g. AF=123; becomes HC_AF=123;), and
#   2. a header file where the info and format fields are
#   updated to match, with prepended caller labels (e.g.
#   INFO=<ID=HC_AF>).
# Note that the GT field is duplicated.  The first entry is left
# as "GT"; the second is prepended with the caller label ("HC_GT").
# This allows bcftools merge to function as intended.


set -euo pipefail

inFile=$1
caller=$2

inFileName="${inFile##*/}"
outFile="prepended.${inFileName}"

# outFile="${saveDir}/prepended.${inFileName}"

cut -f8 "${inFile}" | awk -v var="${caller}" 'BEGIN{FS=OFS=";"} {for(i=1;i<=NF;i++) $i=var"_"$i; print $0}' > "info.${inFileName}.temp"
cut -f9 "${inFile}" | awk -v var="${caller}" 'BEGIN{FS=OFS=":"} {for(i=1;i<=NF;i++) $i=var"_"$i; print "GT:"$0}' > "format.${inFileName}.temp"
cut -f1-7 "${inFile}" > "left.${inFileName}.temp"
cut -f10- "${inFile}" | awk 'BEGIN{FS=OFS="\t"} {for(i=1;i<=NF;i++) {split($i,a,":"); $i=a[1]":"$i} print $0}' > "right.${inFileName}.temp"
    # need to duplicate GT and preserve one as plain "GT" to allow bcftools merge to work as expected
paste "left.${inFileName}.temp" "info.${inFileName}.temp" "format.${inFileName}.temp" "right.${inFileName}.temp" | awk -v var="${caller}" 'BEGIN{FS=OFS="\t"} {$8 = $8";"var"_QUAL="$6; print $0}' 

# paste "left.${inFileName}.temp" "info.${inFileName}.temp" "format.${inFileName}.temp" "right.${inFileName}.temp" | awk -v var="${caller}" 'BEGIN{FS=OFS="\t"} {$8 = $8";"var"_QUAL="$6; print $0}' > "${outFile}"
# rm "left.${inFileName}.temp" "info.${inFileName}.temp" "format.${inFileName}.temp" "right.${inFileName}.temp"
# cat $outFile

# test : add copying qual into INFO field as caller_qual=
