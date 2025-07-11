#!/bin/bash

# Usage: ./prepend_INFO_labels_variant.sh <caller_label>
#
# Takes input from stdin and prepends the caller label to INFO and FORMAT fields.
# It also duplicates the GT field, prepending it with the caller label (e.g. "HC_GT").

set -euo pipefail

# Get the caller label from the second argument
caller=$1

# Process input from stdin using awk
awk -v caller="${caller}" 'BEGIN{FS=OFS="\t"} {

    # Check if the line starts with "#" (header line)
    if ($0 ~ /^#/) {
        header = $0;

        # Only replace ##INFO or ##FORMAT if the line starts with those
        if ($0 ~ /^##INFO/) {
            gsub(/^##INFO=<ID=/, "##INFO=<ID=" caller "_", header);
        }

        if ($0 ~ /^##FORMAT/) {
            gsub(/^##FORMAT=<ID=/, "##FORMAT=<ID=" caller "_", header);
        }
        
        if (header ~ /^##FORMAT=<ID='${caller}_GT/') {
            header = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Concordant Genotype\">\n##INFO=<ID=" caller "_QUAL,Number=1,Type=String,Description=\"Record of original QUAL value\">" "\n" header;
        }

        # Print the modified header line
        print header
        next;  # Skip further processing for header lines
    }

    # Process fields 1-7 as left part
    left = $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7;
    
    # Process field 8 with the caller prefix
    split($8, info, ";");
    for (i = 1; i <= length(info); i++) info[i] = caller"_"info[i];
    info_str = join(info, ";");
    info_str = info_str";"caller"_QUAL=" $6;
    
    # Process field 9 with the caller prefix
    split($9, format, ":");
    for (i = 1; i <= length(format); i++) format[i] = caller"_"format[i];
    format_str = "GT:"join(format, ":");
    
    # Process fields 10 onwards (right part), and modify GT value if necessary
    right = "";
    for (i = 10; i <= NF; i++) {
        split($i, a, ":");
        $i = a[1]":"$i;
        right = right (right == "" ? $i : "\t" $i);
    }
    
    # Add modified QUAL value to the 8th column
    
    # Print the final combined line
    print left"\t"info_str"\t"format_str"\t"right;
} 
# Function to join an array with a separator
function join(arr, sep) {
    str = arr[1];
    for (i = 2; i in arr; i++) str = str sep arr[i];
    return str;
}' 
