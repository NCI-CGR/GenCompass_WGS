#!/usr/bin/env python3

import pysam
import click
import os
import collections
import pandas as pd

ALLOWED_CHROMOSOMES = {'chr' + str(i) for i in range(1, 23)}
ALLOWED_CHROMOSOMES.update(['chrX', 'chrY'])

@click.command()
@click.option('--original_bam_path', type=click.Path(exists=True))
@click.option('-odir','--output_dir', type=click.Path())
def extract_nonprimary_reads(original_bam_path, output_dir):
    sample_id = os.path.basename(original_bam_path).rstrip('.bam')

    processed_recovered_qnames = set()
    output_qname_rg_path = os.path.join(output_dir, f'{sample_id}.qname_rg_map.tsv')

    
    output_bam_path = os.path.join(output_dir, f'{sample_id}.non_primary.bam')

    non_primary_counts = collections.defaultdict(int)
    # Open the second BAM file for reading
    with pysam.AlignmentFile(original_bam_path, "rb") as bam_in:
        # Open the output BAM file for writing
        with pysam.AlignmentFile(output_bam_path, "wb", header=bam_in.header) as bam_out, open(output_qname_rg_path, 'w') as qname_rg_out:
            qname_rg_out.write(f'QNAME\tRG\n')
            for read in bam_in.fetch(until_eof=True):
                # Process each read
                # Check if the read or its mate is mapped to a primary chromosome
                read_chrom = read.reference_name
                mate_chrom = read.next_reference_name
                if read.is_unmapped or read_chrom not in ALLOWED_CHROMOSOMES or mate_chrom not in ALLOWED_CHROMOSOMES:
                    bam_out.write(read)
                    non_primary_counts[read_chrom] = non_primary_counts[read_chrom] + 1
                    if read.query_name not in processed_recovered_qnames:
                        processed_recovered_qnames.add(read.query_name)
                        qname_rg_out.write(f'{read.query_name}\t{read.get_tag("RG")}\n')

                    # qname_rg_map[read.query_name] = read.get_tag('RG')
    
    output_qname_rg_path = os.path.join(output_dir, f'{sample_id}.qname_rg_map.tsv')
    # qname_rg_df = pd.DataFrame(list(qname_rg_map.items()), columns=['QNAME', 'RG'])

    # qname_rg_df.to_csv(output_qname_rg_path, sep='\t', index=False)


    for key, value in non_primary_counts.items():
        print(f'{key}\t{value}')

if __name__ == "__main__":
    extract_nonprimary_reads()