#!/usr/bin/env python3
import pandas as pd
import os
import argparse
import json

def add_read_group(row, platform='ILLUMINA'):
    date_str = pd.to_datetime(row["SEQDATE"]).strftime("%Y%m%d")
    read_group = f"@RG\\tID:{row['FLOWCELL']}.{row['LANE']}\\tLB:{row['CGF ID']}_{row['INDEX']}\\tPL:{platform}\\tSM:{row['GROUP']}_{row['ANALYSIS ID']}\\tPU:{row['FLOWCELL']}{date_str}.{row['LANE']}.{row['INDEX']}\\tCN:CGR"
    return read_group

def add_fastq_substring(row):
    fastq_substring = f"{row['CGF ID']}_{row['INDEX']}_L00{row['LANE']}_R1" 
    return fastq_substring

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', help='project manfiest file. Contains "ANALYSIS ID", "FLOWCELL", "INDEX", "" columns')
    parser.add_argument('--fastq_files', nargs='*', help='Fastq filenames being converted to BAM files')
    parser.add_argument('--analysis_id', help='CGR Analysis ID')
    parser.add_argument('-odir', '--output_directory', default="./")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    
    if len(args.fastq_files) > 1: 
        fastq_files = pd.DataFrame(args.fastq_files, columns=['Fastq File'])
    else:
        fastq_files = pd.read_csv(args.fastq_files[0], header=None, names=['Fastq File'])
    
    manifest = pd.read_csv(args.manifest)
    manifest['RG'] = manifest.apply(add_read_group, axis=1)
    manifest['Fastq_Substring'] = manifest.apply(add_fastq_substring, axis=1)
    
    match = manifest.loc[manifest["ANALYSIS ID"] == args.analysis_id]
    
    if match.empty:
        raise ValueError(f"analysis_id {args.analysis_id} not found")
    
    # Create map of destination -> source for fastq files
    fastq_map = {}
    fq2bam_entries = []
    
    # Process each matching row
    for _, manifest_row in match.iterrows():
        read_group = manifest_row['RG']
        fastq_subtr = manifest_row['Fastq_Substring']
        flowcell = manifest_row['FLOWCELL']
        lane = manifest_row['LANE']
        
        # Create unique directory for this flowcell.lane combination
        unique_dir = f"{flowcell}.{lane}"
        
        filtered = fastq_files[
            fastq_files["Fastq File"].str.contains(fastq_subtr, na=False) &
            fastq_files["Fastq File"].str.contains(flowcell, na=False) &
            fastq_files["Fastq File"].str.contains(f'/L{lane}/', na=False)
        ]
        
        
        # Add filtered fastq files to the map with unique paths
        for _, fastq_row in filtered.iterrows():
            r1_source = fastq_row['Fastq File']
            r1_basename = os.path.basename(r1_source)
            r2_basename = r1_basename.replace('_R1_001_HQ_paired.fastq.gz', '_R2_001_HQ_paired.fastq.gz')
            r2_source = r1_source.replace('_R1_001_HQ_paired.fastq.gz', '_R2_001_HQ_paired.fastq.gz')
            
            # Create mapped paths
            r1_dest = f"{unique_dir}/{r1_basename}"
            r2_dest = f"{unique_dir}/{r2_basename}"
            
            # Add to map
            fastq_map[r1_dest] = r1_source
            fastq_map[r2_dest] = r2_source
            
            # Add to fq2bam list with new paths
            fq2bam_entries.append({
                'r1': r1_dest,
                'r2': r2_dest,
                'rg': read_group
            })
    
    # Write fq2bam list with mapped paths
    with open(os.path.join(args.output_directory, f'{args.analysis_id}.fq2bam_list.txt'), 'w') as ofile:
        for entry in fq2bam_entries:
            ofile.write(f"{entry['r1']}\t{entry['r2']}\t\"{entry['rg']}\"\n")
    
    # Write fastq map as JSON
    with open(os.path.join(args.output_directory, f'{args.analysis_id}_fastq_map.json'), 'w') as ofile:
        json.dump(fastq_map, ofile, indent=2)
    
    # Write original fastq files list for reference
    with open(os.path.join(args.output_directory, f'{args.analysis_id}_fastq_files.txt'), 'w') as ofile:
        for source in set(fastq_map.values()):
            ofile.write(f'{source}\n')

if __name__=="__main__":
    main()