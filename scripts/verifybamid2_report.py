#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import datetime
from helpers import build_filelist


def save_verifybamid_selfsm_report(report, odir, project):
    column_descriptions = pd.DataFrame(
        [
            ["#SEQ_ID","Sample ID"],
            ["#SNPS", "# of SNPs passing the criteria from the VCF file"],
            ["#READS", "Total # of reads loaded from the BAM file"],
            ["AVG_DP", "Average sequencing depth at the sites in the VCF file"],
            ["FREEMIX", "Estimate of contamination (0-1 scale) based on sequencing data only (not sequence+array)"]
        ],
        columns=['Column', 'Description']
    )

    current_date = datetime.datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"
   
    filename = os.path.join(odir, f'{project}.verifybamid2.selfSM-{date_suffix}.xlsx')
    with pd.ExcelWriter(filename) as writer:
        report.to_excel(writer, sheet_name='VerifyBamID2.selfSM Report', index=False)
        column_descriptions.to_excel(writer, sheet_name="Column Descriptions", index=False)

        workbook = writer.book
        worksheet = workbook.get_worksheet_by_name('VerifyBamID2.selfSM Report')
        thousands_format = workbook.add_format({'num_format': '#,###'})
        
        thousands_columns = ["#SNPS", "#READS",]
        for col in thousands_columns:
            index_no = report.columns.get_loc(col)
            worksheet.set_column(index_no, index_no, None, thousands_format)

def read_selfsm(filename):
    selfsm = pd.read_csv(filename, sep='\t', 
        usecols=[
        "#SEQ_ID",
        "#SNPS",
        "#READS",
        "AVG_DP",
        "FREEMIX"
        ])
    return selfsm

def create_verifybamid_selfsm_report(selfsm_files):
    verifybamid_report = [read_selfsm(filename) for filename in selfsm_files]
    verifybamid_report = pd.concat(verifybamid_report, ignore_index=True)
    return verifybamid_report


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--selfsm', nargs='*', required=False, help='SelfSM files generated from verifybamid2')
    parser.add_argument('-sl', '--selfsm_list', required=False, help='List of SelfSM files generated from verifybamid2')
    parser.add_argument('-odir', '--output_directory', dest='odir', default='./')
    parser.add_argument('-p', '--project', help='project name')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    selfsm_list = build_filelist(args.selfsm, args.selfsm_list)
    verifybamid_selfsm_report = create_verifybamid_selfsm_report(selfsm_list)
    
    save_verifybamid_selfsm_report(verifybamid_selfsm_report, args.odir, args.project)


if __name__=="__main__":
    main()
    
    
