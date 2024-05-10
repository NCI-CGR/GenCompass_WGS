#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import csv
import helpers


class Fq2BamInputBuilder:
    def __init__(self, fastq_locations: pd.DataFrame, manifest: helpers.Manifest, sample_id: str, platform:str) -> None:
        self.fastq_locations = fastq_locations.copy()

        self.filename_column = self.fastq_locations.columns[-1]
        self.fastq_locations['Basename'] = self.fastq_locations[self.filename_column].apply(lambda fpath: os.path.basename(fpath))
        self.sample_id = sample_id
        self.manifest = manifest
        self.platform = platform
        
        self.fastq_locations['Sample Run ID'] = self.fastq_locations[self.filename_column].apply(helpers.extract_sample_run_id)
        print(self.fastq_locations)

        self.filter_sample_locations()
        self.sample_paired_end_table = self.build_sample_paired_end_table()
        self._add_read_group_to_paired_end()


    def build_sample_paired_end_table(self):
        self.fastq_locations['Sample Lane ID'] = self.fastq_locations['Basename'].apply(helpers.extract_sample_lane_id)
        self.fastq_locations['Paired End'] = self.fastq_locations['Basename'].apply(helpers.extract_sample_paired_end)
        sample_paired_end_table =self.fastq_locations.pivot(
            index=['Sample Run ID','Sample Lane ID' ],
            columns='Paired End', 
            values='Basename')
        sample_paired_end_table = sample_paired_end_table[["R1", "R2"]]
        sample_paired_end_table.reset_index(inplace=True)
        return sample_paired_end_table
        

    def filter_sample_locations(self):
        sample_run_ids = self.manifest.get_sample_run_ids(self.sample_id)
        self.fastq_locations = self.fastq_locations[self.fastq_locations['Sample Run ID'].isin(sample_run_ids)].copy()

    # def _extract_sample_run_id(self, fastq_path: str) -> str:
    #     sample_run_id = os.path.basename(fastq_path).split('_')[0]
    #     return sample_run_id

    # def _extract_sample_lane_id(self, fastq_path: str) -> str:
    #     """Extracts the sample lane ID from the fastq filename.
    #     Lane ID is all parts of the filename before the paired end

    #     Example:
    #         input: path/to/files/I3-98765_S23_L001_R1_001.fastq.gz
    #         output: I3-98765_S23_L001

    #         input: path/to/files/I3-98765_S23_L001.R1.001.fastq.gz
    #         output: I3-98765_S23_L001

    #         input: path/to/files/I3-98765_S23_L001.R1_001.fastq.gz
    #         output: I3-98765_S23_L001

    #         input: path/to/files/I3-98765_S23_L001_R1.fastq.gz
    #         output: I3-98765_S23_L001

    #     Args:
    #         fastq_path (str): path to the fastq file

    #     Returns:
    #         str: Lane ID
    #     """
    #     fastq_filename = os.path.basename(fastq_path)
    #     lane_id_regex = re.compile(r'^.*(?=((\.|\_)R\d(\.|\_)))')
    #     lane_id = lane_id_regex.search(fastq_filename)
    #     return 'N/A' if lane_id is None else lane_id.group(0)


    # def _extract_sample_paired_end(self, fastq_path: str) -> str:
    #     """Extracts the paired end of the fastq file.
    #     The paired end is epected to be R1 or R2, and separated
    #     from the rest of the filename by either a "." or "_"
    #     Example:
    #         input: path/to/files/I3-98765_S23_L001_R1_001.fastq.gz
    #         output: R1

    #         input: path/to/files/I3-98765_S23_L001.R1.001.fastq.gz
    #         output: R1

    #         input: path/to/files/I3-98765_S23_L001.R1_001.fastq.gz
    #         output: R1

    #         input: path/to/files/I3-98765_S23_L001_R1.fastq.gz
    #         output: R1


    #     Args:
    #         fastq_path (str): path of the fastq file

    #     Returns:
    #         str: paired end string
    #     """
    #     fastq_filename = os.path.basename(fastq_path)
    #     paired_end_regex= re.compile(r'(\.|\_)R\d(\.|\_)')
    #     paired_end = paired_end_regex.search(fastq_filename)
    #     if paired_end is not None:
    #         paired_end = paired_end.group(0)
    #         paired_end = paired_end.strip('.').strip('_')
    #         return paired_end
    #     else:
    #         return "N/A"

    def _add_read_group_to_paired_end(self):
        def get_read_group(row):
            sample, index, lane = row['Sample Lane ID'].split('_')
            return f"@RG\\tID:{row['Sample Lane ID']}\\tSM:{row[self.manifest.sample_id_column]}\\tPL:{self.platform}\\tLB:{row['Sample Lane ID']}\\tPU:{row[self.manifest.flowcell_column]}.{lane}.{index}"
        self.sample_paired_end_table = self.sample_paired_end_table.merge(self.manifest.manifest, left_on='Sample Run ID', right_on=self.manifest.sample_run_id_column)
        # self.sample_paired_end_table['Sample ID'] = self.sample_paired_end_table['Sample Run ID'].apply(lambda run: self.manifest.get_sample_id(run))
        self.sample_paired_end_table['Read Group'] = self.sample_paired_end_table.apply(get_read_group, axis=1)
        self.sample_paired_end_table['Read Group'] = self.sample_paired_end_table['Read Group'].apply(lambda rg: '"' + rg + '"')
        self.sample_paired_end_table = self.sample_paired_end_table[['R1', 'R2', 'Read Group']].copy()

    def save_read_group_file(self, odir):
        ofilename = os.path.join(odir, f'{self.sample_id}.fq2bam_list.txt')
        self.sample_paired_end_table.to_csv(ofilename, sep=' ', index=False, header=False, quoting=csv.QUOTE_NONE)     

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', help='project manfiest file. Contains "Sample Run ID", "Sample ID", "Flowcell" columns')
    parser.add_argument('--fastq_files', nargs='*', help='Fastq filenames being converted to BAM files')
    parser.add_argument('--platform', action='store', default='ILLUMINA')
    parser.add_argument('-odir', '--output_directory', default="./")
    parser.add_argument('--sample', default="sample")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    # fastq_files = read_table(args.fastq_files, header=None)
    if len(args.fastq_files) > 1: 
        fastq_files = pd.DataFrame(args.fastq_files, columns=['Fastq File'])
    else:
        fastq_files = pd.read_csv(args.fastq_files[0], header=None, names=['Fastq File'])
    
    manifest = helpers.Manifest(args.manifest)
    fq2bam_builder = Fq2BamInputBuilder(fastq_files, manifest, args.sample, args.platform)
    fq2bam_builder.save_read_group_file(args.output_directory)
    fq2bam_builder.fastq_locations[['Fastq File']].to_csv(os.path.join(args.output_directory,f'{args.sample}_fastq_files.txt'), index=False, header=False)

if __name__=="__main__":
    main()