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

headerFile=$1
caller=$2
saveDir=$3

headerName="${headerFile##*/}"
outHeader="${saveDir}/prepended.${headerName}"

sed "s/##INFO=<ID=/##INFO=<ID=${caller}_/" "${headerFile}" | sed "s/##FORMAT=<ID=/##FORMAT=<ID=${caller}_/" | sed "/##FORMAT=<ID=${caller}_GT/ i\##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Concordant Genotype\">\n##INFO=<ID=${caller}_QUAL,Number=1,Type=String,Description=\"Record of original QUAL value\">" > "${outHeader}"

# test : add copying qual into INFO field as caller_qual=
