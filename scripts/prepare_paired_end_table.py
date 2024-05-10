#!/usr/bin/env python3

"""
Creates a table of fastq locations of a given sample formatted for WDL workflow
Each Row is a lane of the sample. There is no header included for the table.
The 

Input:
| Argument              | Description                                                       |
|-----------------------|-------------------------------------------------------------------|
| -i, --fastq_locations | Table of fastq locations (.csv, .txt, .tsv, .xls or .xlsx format) |
| -id, --sample_id      | Sample ID                                                         |
| -odir                 | Output Directory                                                  |

Output:
SAMPLE_ID.fastq_locations.tsv : Outputs tab separated table without headers

-----------------
Example

Input:
sample_id: SC100000

Manifest
| Sample ID | Sample Run ID | Other Columns... |
| ----------| --------------|------------------|
| SC100000  | I3-12345      | ....             |
| SC100000  | I3-97534      | ....             |

Fastq Locations 

path/to/files/I3-98765_S23_L001_R1_001.fastq.gz
path/to/files/I3-98765_S23_L001_R2_001.fastq.gz
path/to/files/I3-98765_S23_L002_R1_001.fastq.gz
path/to/files/I3-98765_S23_L002_R2_001.fastq.gz
---
or
---
| Sample ID | Fastq Path |
| --------- | --------------|
| SC100000  | path/to/files/I3-12345_S23_L001_R1_001.fastq.gz |
| SC100000  | path/to/files/I3-12345_S23_L001_R1_001.fastq.gz |
| SC100000  | path/to/files/I3-12345_S23_L001_R2_001.fastq.gz |
| SC100000  | path/to/files/I3-12345_S23_L002_R1_001.fastq.gz |
| SC100000  | path/to/files/I3-12345_S23_L002_R2_001.fastq.gz |

Example Output (headers not included in output):
| Sample ID      |  Lane ID          | R1 file location                                | R2 file location                                |
|----------------|-------------------|-------------------------------------------------|-------------------------------------------------|
| SC100000       | I3-98765_S23_L001 | path/to/files/I3-98765_S23_L001_R1_001.fastq.gz | path/to/files/I3-98765_S23_L001_R2_001.fastq.gz |
| SC100000       | I3-98765_S23_L002 | path/to/files/I3-98765_S23_L002_R1_001.fastq.gz | path/to/files/I3-98765_S23_L002_R2_001.fastq.gz |

"""

import pandas as pd
import os
import re
import argparse
import helpers


class SampleFileBuilder:
    def __init__(self, fastq_locations_file: str, manifest: helpers.Manifest, sample_id: str, omics_fastq_map_file: str = None, aws_account_id: str = None) -> None:

        self.sample_id = sample_id
        self.manifest = manifest
        self.aws_account_id=aws_account_id
        
        if omics_fastq_map_file is None:
            self.fastq_locations = helpers.read_table(fastq_locations_file, header=None)
            self.filename_column = self.fastq_locations.columns[-1]
            self.fastq_locations['Sample Run ID'] = self.fastq_locations[self.filename_column].apply(helpers.extract_sample_run_id)
            self.filter_sample_locations()
            self.fastq_locations['Sample ID'] = sample_id
            self.sample_paired_end_table = self.build_sample_paired_end_table()
        else:
            self.read_omics_fastq_map(omics_fastq_map_file)
            self.sample_paired_end_table = self.build_omics_sample_paired_end_table()

    def build_omics_sample_paired_end_table(self):
        sample_run_ids = self.manifest.get_sample_run_ids(self.sample_id)
        self.omics_fastq_map =  self.omics_fastq_map[self.omics_fastq_map['Sample Run ID'].isin(sample_run_ids)].copy()
        self.omics_fastq_map['Sample ID'] = self.sample_id

        sample_paired_end_table =self.omics_fastq_map.pivot(
            index=['Sample ID', 'Sample Lane ID',], 
            columns='Paired End', 
            values='Omics Storage Location')
        sample_paired_end_table = sample_paired_end_table[["R1", "R2"]]
        sample_paired_end_table.reset_index(inplace=True)
        return sample_paired_end_table


    def build_sample_paired_end_table(self):
        self.fastq_locations['Sample Lane ID'] = self.fastq_locations[self.filename_column].apply(helpers.extract_sample_lane_id)
        self.fastq_locations['Paired End'] = self.fastq_locations[self.filename_column].apply(helpers.extract_sample_paired_end)
        sample_paired_end_table =self.fastq_locations.pivot(
            index=['Sample ID', 'Sample Lane ID',], 
            columns='Paired End', 
            values=self.filename_column)
        sample_paired_end_table = sample_paired_end_table[["R1", "R2"]]
        sample_paired_end_table.reset_index(inplace=True)
        return sample_paired_end_table
        
    def read_omics_fastq_map(self, omics_fastq_map_file):
        self.omics_fastq_map=pd.read_csv(omics_fastq_map_file, sep='\t', header=None, names=['omics_id', 'fastq'])
        
        self.omics_fastq_map['Sample Run ID'] = self.omics_fastq_map['fastq'].apply(lambda x: x.split('_')[0])
        self.omics_fastq_map['Paired End'] = self.omics_fastq_map['fastq'].apply(lambda x: x.split('_')[3])
        self.omics_fastq_map['Sample Lane ID'] = self.omics_fastq_map['fastq'].apply(lambda x: '_'.join(x.split('_')[:3]))
        self.omics_fastq_map['Sequence Store ID'] = self.omics_fastq_map['omics_id'].apply(lambda x: x.split('_')[0])
        self.omics_fastq_map['Read Set ID'] = self.omics_fastq_map['omics_id'].apply(lambda x: x.split('_')[1])

        def get_omics_storage_location(row):
            storage_map = {'R1' : 'source1', 'R2': 'source2'}
            omics_storage_location = f"omics://{self.aws_account_id}.storage.us-east-1.amazonaws.com/{row['Sequence Store ID']}/readSet/{row['Read Set ID']}/{storage_map[row['Paired End']]}"
            return omics_storage_location
        self.omics_fastq_map['Omics Storage Location'] = self.omics_fastq_map.apply(get_omics_storage_location, axis=1)

    def filter_sample_locations(self):
        sample_run_ids = self.manifest.get_sample_run_ids(self.sample_id)
        self.fastq_locations = self.fastq_locations[self.fastq_locations['Sample Run ID'].isin(sample_run_ids)].copy()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fastq_locations', action='store', dest='fastq_locations', help='table that has the fastq locations', required=False)
    parser.add_argument('-m', '--manifest', action='store', dest='manifest', help='Manifest table')
    parser.add_argument('-id','--sample_id', action='store', dest='sample_id')
    parser.add_argument('-omap', '--omics_fastq_map_file', action='store', dest='omics_fastq_map_file', required=False)
    parser.add_argument('-aws_acct', '--aws_account_id',  action='store', dest='aws_account_id', required=False)
    parser.add_argument('-odir', action='store', dest='odir', default='./')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    os.makedirs(args.odir, exist_ok = True) 

    manifest = helpers.Manifest(args.manifest)
    
    sample = SampleFileBuilder(args.fastq_locations, manifest, sample_id=args.sample_id, omics_fastq_map_file=args.omics_fastq_map_file, aws_account_id=args.aws_account_id)
    filepath = os.path.join(args.odir, f'{args.sample_id}.fastq_locations.tsv')
    sample.sample_paired_end_table.to_csv(filepath, sep='\t', index=False, header=False)



if __name__=="__main__":
    main()
    
