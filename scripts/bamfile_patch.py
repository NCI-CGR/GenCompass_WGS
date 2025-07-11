#!/usr/bin/env python3

import pysam
import click
import os
from collections import defaultdict
import time
import logging
import subprocess
import pandas as pd

logger = logging.getLogger(__name__)

ALLOWED_CHROMOSOMES = {'chr' + str(i) for i in range(1, 23)}
ALLOWED_CHROMOSOMES.update(['chrX', 'chrY'])

@click.group()
def cli():
    pass

def replace_read_group(tags, read_group):
    # tags.append(('PG', 'BamPatch'))
    return [('RG', read_group) if tag == 'RG' else (tag, value) for tag, value in tags]

def build_patch_qnames(bam_file):
    logger.info('Building patch qnames set')
    qnames = set()
    start = time.time()
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            qname = read.query_name
            qnames.add(qname)
    end = time.time()
    logger.info(f'Time to build qname set \t: {end - start:.2f} seconds')
    return qnames
    
@cli.command()
@click.option('--non_primary_bam_path', type=click.Path(exists=True))
@click.option('--original_bam_path', type=click.Path(exists=True))
@click.option('--qname_rg_map', type=click.Path(exists=True))
@click.option('--sample_id', type=click.STRING)
@click.option('-odir','--output_dir', type=click.Path(), default='./')
def patch_bam(non_primary_bam_path,original_bam_path, qname_rg_map, sample_id, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    # sample_id = os.path.basename(original_bam_path).rstrip('.bam')
    logging.basicConfig(filename=f'{output_dir}/{sample_id}.bam_patch.log', level=logging.INFO)
    start = time.time()

    qname_rg_df = pd.read_csv(qname_rg_map, sep='\t')
    original_read_groups = pd.Series(qname_rg_df['RG'].values, index=qname_rg_df['QNAME']).to_dict()

    end = time.time()
    lap = end
    
    mark_duplicates_pg_entry = {
    'ID': 'MarkDuplicates',
    'PN': 'MarkDuplicates',
    'PP': 'pbrun fq2bam',
    'VN': '4.3.0.0',
    'CL': 'pbrun fq2bam'
    }

    bamfile_patch_pg_entry = {
    'ID': 'BamPatch',
    'PN': 'BamPatch',
    'PP': 'MarkDuplicates',
    'VN': '1.2',
    'CL': f'bamfile_patch.py'
    }

    original_bam = pysam.AlignmentFile(original_bam_path, "rb", )
    base_header = original_bam.header.to_dict()
    original_bam.close()

    pg_entries = base_header.get('PG', [])
    base_rg = base_header.get('RG', [])
    pg_entries.append(mark_duplicates_pg_entry)
    # pg_entries.append(bamfile_patch_pg_entry)
    base_header['PG'] = pg_entries    
    original_bam.close()

    patch_out_path = os.path.join(output_dir, f'{sample_id}.patch_rg_update.bam')
    patch_bam = pysam.AlignmentFile(non_primary_bam_path, "rb")
    patch_header = patch_bam.header.to_dict()
    pg_entries = patch_header.get('PG', [])
    pg_entries.append(mark_duplicates_pg_entry)
    pg_entries.append(bamfile_patch_pg_entry)
    patch_header['RG'] = base_rg
    patch_header['PG'] = pg_entries
    with pysam.AlignmentFile(patch_out_path, "wb", header=patch_header) as patch_out:
        logger.info(f'\n{"*"*75}\nUpdating patch with original RG')
        for read in patch_bam.fetch(until_eof=True):
            qname = read.query_name
            read.set_tags(replace_read_group(read.get_tags(), original_read_groups[qname]))
            patch_out.write(read)
    patch_bam.close()
    end = time.time()
    logger.info(f'Time to process patch RG update: {end - lap: .2f} seconds. Total time {(end - start) / 60: .2f} minutes' )

@cli.command()
@click.option('--non_primary_bam_path', type=click.Path(exists=True))
@click.option('--sample_id', type=click.STRING)
@click.option('-odir','--output_dir', type=click.Path(), default='./')
def qname(non_primary_bam_path, sample_id, output_dir):
    patch_qnames = build_patch_qnames(non_primary_bam_path)
    with open(f'{output_dir}/{sample_id}_qnames.txt', 'w') as ofile:
        for qname in patch_qnames:
            ofile.write(f"{qname}\n")

if __name__ == "__main__":
    cli()
