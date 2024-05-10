#!/usr/bin/env python3

import pandas as pd
import os
import re 

class Manifest:
    def __init__(self, manifest_file: str) -> None:
        self.manifest_file = manifest_file
        self.manifest = read_table(self.manifest_file)
        self.sample_id_column = self.get_column_match('Sample ID')
        self.sample_run_id_column = self.get_column_match('Sample Run ID')
        if self.sample_run_id_column is None:
            self.sample_run_id_column = self.sample_id_column

        self.flowcell_column = self.get_column_match('Flowcell')   
        
    
    def get_column_match(self, column: str):
        standard_column = standardize_string(column)
        for col in self.manifest.columns:
            std_manifest_column = standardize_string(col)
            if std_manifest_column == standard_column:
                return col
        return None
    
    def get_sample_run_ids(self, sample_id):
        sample: pd.DataFrame = self.manifest[self.manifest[self.sample_id_column] == sample_id]
        sample_run_ids = sample[self.sample_run_id_column].unique()
        return list(sample_run_ids)
    
    def get_sample_run_id_map(self):
        sample_run_map = pd.Series(self.manifest[self.sample_id_column].values,index=self.manifest[self.sample_run_id_column]).to_dict()
        return sample_run_map
        
        

def read_table(table_filename: str, sheet_name = 0, header=0) -> pd.DataFrame:
    extension = os.path.splitext(table_filename)[1]
    extension = extension.lower()
    if extension in (".tsv", ".txt"):
        return pd.read_csv(table_filename, sep = '\t', header=header)
    elif extension in (".csv"):
        return pd.read_csv(table_filename, header=header)
    elif extension in (".xlsx", ".xls"):
        return pd.read_excel(table_filename, sheet_name=sheet_name, header=header)
    else:
        raise ValueError("Unknown file type. File type must be one of [.csv, .txt, .tsv, .xls, .xlsx]")
    

def standardize_string(string: str) -> str:
        standard_string = string.replace(' ', '')
        standard_string = standard_string.replace('_', '')
        standard_string = standard_string.lower()
        return standard_string



def extract_sample_run_id(fastq_path: str) -> str:
    sample_run_id = os.path.basename(fastq_path).split('_')[0]
    return sample_run_id

def extract_sample_lane_id(fastq_path: str) -> str:
    """Extracts the sample lane ID from the fastq filename.
    Lane ID is all parts of the filename before the paired end

    Example:
        input: path/to/files/I3-98765_S23_L001_R1_001.fastq.gz
        output: I3-98765_S23_L001

        input: path/to/files/I3-98765_S23_L001.R1.001.fastq.gz
        output: I3-98765_S23_L001

        input: path/to/files/I3-98765_S23_L001.R1_001.fastq.gz
        output: I3-98765_S23_L001

        input: path/to/files/I3-98765_S23_L001_R1.fastq.gz
        output: I3-98765_S23_L001

    Args:
        fastq_path (str): path to the fastq file

    Returns:
        str: Lane ID
    """
    fastq_filename = os.path.basename(fastq_path)
    lane_id_regex = re.compile(r'^.*(?=((\.|\_)R\d(\.|\_)))')
    lane_id = lane_id_regex.search(fastq_filename)
    return 'N/A' if lane_id is None else lane_id.group(0)


def extract_sample_paired_end(fastq_path: str) -> str:
    """Extracts the paired end of the fastq file.
    The paired end is epected to be R1 or R2, and separated
    from the rest of the filename by either a "." or "_"
    Example:
        input: path/to/files/I3-98765_S23_L001_R1_001.fastq.gz
        output: R1

        input: path/to/files/I3-98765_S23_L001.R1.001.fastq.gz
        output: R1

        input: path/to/files/I3-98765_S23_L001.R1_001.fastq.gz
        output: R1

        input: path/to/files/I3-98765_S23_L001_R1.fastq.gz
        output: R1


    Args:
        fastq_path (str): path of the fastq file

    Returns:
        str: paired end string
    """
    fastq_filename = os.path.basename(fastq_path)
    paired_end_regex= re.compile(r'(\.|\_)R\d(\.|\_)')
    paired_end = paired_end_regex.search(fastq_filename)
    if paired_end is not None:
        paired_end = paired_end.group(0)
        paired_end = paired_end.strip('.').strip('_')
        return paired_end
    else:
        return "N/A"


def build_filelist(files, files_fof):
    files = [] if files is None else files
    with open(files_fof) as f:
        files_from_fof =  f.read().splitlines()
    files = files + files_from_fof
    return files