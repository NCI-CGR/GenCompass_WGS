#!/usr/bin/env python3

import argparse
import pandas as pd
import os
from helpers import build_filelist

def extract_bammmetric_data(bammetrics_filename: str):
    
    bammetrics = pd.read_csv(bammetrics_filename,  sep='\t', skiprows=1, nrows=1)
    sample_name = os.path.basename(bammetrics_filename).strip('.bammetrics.txt')
    bammetrics.insert(0,column='Sample', value=sample_name)
    return bammetrics


def create_bammetrics_report(bammetrics_files: list):
    bammetrics_report = [extract_bammmetric_data(filename) for filename in bammetrics_files]
    bammetrics_report = pd.concat(bammetrics_report)
    return bammetrics_report


def save_report(report, output_directory, project):
    output_filename = os.path.join(output_directory, f'{project}.bammetrics_report.tsv')
    report.to_csv(output_filename,  index=False, sep='\t')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--bammetrics', nargs='*',  help='bammetrics file(s) generated by parabricks', required=False)
    parser.add_argument('-l', '--bammetrics_list', required = False, help='List of bammetrics files generated by parabricks')
    parser.add_argument('-odir',default = './', help='Output directory to save bammetrics report')
    parser.add_argument('--project', help='Name of project. Used for saving file project.bammetrics_report.xlsx')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    bammetrics_list = build_filelist(args.bammetrics, args.bammetrics_list)

    report = create_bammetrics_report(bammetrics_list)
    save_report(report, args.odir, args.project)

if __name__=='__main__':
    main()
