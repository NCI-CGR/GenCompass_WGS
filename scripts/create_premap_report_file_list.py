#!/usr/bin/env python3

import pandas as pd
import helpers
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--sample_list')
parser.add_argument('--manifest')
parser.add_argument('--sample_fastq_files')
parser.add_argument('--workflow_results_dir')


args = parser.parse_args()

manifest = helpers.Manifest(args.manifest)

with open(args.sample_list) as f:
    samples =  f.read().splitlines()

manifest.manifest = manifest.manifest[manifest.manifest[manifest.sample_id_column].isin(samples)]

fastq_files = helpers.build_filelist(files=None, files_fof=args.sample_fastq_files)
fastq_files = pd.DataFrame(fastq_files, columns=['Path'])
fastq_files['RunID'] = fastq_files['Path'].apply(helpers.extract_sample_run_id)
fastq_files['LaneID'] = fastq_files['Path'].apply(helpers.extract_sample_lane_id)
fastq_files['BaseID'] = fastq_files['Path'].apply(lambda path: os.path.basename(path).rstrip('.fastq.gz'))

with open('fastp_file_list.txt', 'w') as fastp_files:
    for sample in samples:
        run_ids = manifest.get_sample_run_ids(sample)
        sample_lane_ids = fastq_files[fastq_files['RunID'].isin(run_ids)]['LaneID'].unique()
        for lane_id in sample_lane_ids:
            fastp_files.write(os.path.join(args.workflow_results_dir,f'premap_qc/fastp/{sample}/{lane_id}_fastp.json') + '\n')

with open('fastq_screen_files.txt', 'w') as fastq_screen_files:
    for sample in samples:
        fastq_screen_files.write(os.path.join(args.workflow_results_dir, f'premap_qc/fastq_screen/{sample}/{sample}_screen.txt') + '\n')


with open('fastqc_files.txt', 'w') as fastqc_files:
    for sample in samples:
        run_ids = manifest.get_sample_run_ids(sample)
        sample_base_ids = fastq_files[fastq_files['RunID'].isin(run_ids)]['BaseID'].unique()
        print(sample_base_ids)
        for base_id in sample_base_ids:
            fastqc_files.write(os.path.join(args.workflow_results_dir,f'premap_qc/fastqc/{sample}/{base_id}_fastqc.zip') + '\n')


