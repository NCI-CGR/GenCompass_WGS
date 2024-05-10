#!/usr/bin/env python3

import os
import pandas as pd
import argparse
import datetime
from helpers import build_filelist

def extract_fastq_screen(fastq_screen_filename: str):
    def get_num_rows():
        num_rows=0
        with open(fastq_screen_filename) as ifile:
            for line in ifile.readlines():
                if line.startswith('%Hit_no_genomes:'):
                    return num_rows - 3

                else:
                    num_rows = num_rows + 1
        return num_rows
        
    sample = os.path.basename(fastq_screen_filename).strip('_screen.txt')
    nrows = get_num_rows()
    
    fastq_screen_data = pd.read_csv(fastq_screen_filename, skiprows=1, nrows=nrows, sep='\t')
    fastq_screen_data.insert(0, 'Sample', sample)
    return fastq_screen_data




def create_fastq_screen_report(fastq_screen_files: list):
    fastq_data = [extract_fastq_screen(filename) for filename in fastq_screen_files]
    fastq_data = pd.concat(fastq_data, ignore_index=True)
    return fastq_data

def save_fastq_screen_report(data, project, odir):
    column_descriptions = pd.DataFrame(
        [
            ["Sample", "Unique CGR Sample ID"],
            ["Genome", "Run mapping against Human, Yeast, Ecoli, PhiX, Lambda, Vectors, Adapters, Cow, Pig genomes"],
            ["#Reads_processed", "Number of R1 reads downsampled to run FASTQ-Screen (default: 10 million)"],
            ["#Unmapped", "Number of R1 reads that did not map to one the genomes above"],
            ["%Unmapped", "% of R1 reads that did not map to one of the genomes above"],
            ["#One_hit_one_genome", "# R1 reads that mapped uniquely to the specified genome"],
            ["%One_hit_one_genome", "% R1 reads that mapped uniquely to the specified genome"],
            ["#Multiple_hits_one_genome", "# R1 reads that multi-mapped to the specified genome"],
            ["%Multiple_hits_one_genome", "% R1 reads that multi-mapped to the specified genome"],
            ["#One_hit_multiple_genomes", "# R1 reads that mapped uniquely to the specified genome and mapped to at least one other genome"],
            ["%One_hit_multiple_genomes", "% R1 reads that mapped uniquely to the specified genome and mapped to at least one other genome"],
            ["Multiple_hits_multiple_genomes", "# R1 reads that multi-mapped to the specified genome and mapped to at least one other genome"],
            ["%Multiple_hits_multiple_genomes", "% R1 reads that multi-mapped to the specified genome and mapped to at least one other genome"],
        ]
    )

    current_date = datetime.datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"
    output_filename = os.path.join(odir, f'{project}.pre-mapping-QC-B-fastq_screen_report-{date_suffix}.xlsx')
    with pd.ExcelWriter(output_filename) as writer:
        data.to_excel(writer, sheet_name='Fastq Screen Report', index=False)
        column_descriptions.to_excel(writer, sheet_name="Column Descriptions", index=False)
        workbook  = writer.book
        # worksheet = writer.sheets['Fastp_Report']
        worksheet = workbook.get_worksheet_by_name('Fastq Screen Report')
        thousands_unit_columns = [col for col in data.columns if col.startswith('#')]
        thousands_format = workbook.add_format({'num_format': '#,##0'})
        
        for col in thousands_unit_columns:
            index_no = data.columns.get_loc(col)
            worksheet.set_column(index_no, index_no, None, thousands_format)

        
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fastq_screen_files', nargs='*', required=False, help='Fastq screen files to combine into single file')
    parser.add_argument('-l','--fastq_screen_list',  required=False, help='Fastq screen files to combine into single file')
    parser.add_argument('-p', '--project_name', help='Name of project')
    parser.add_argument('-odir','--output_directory', dest = 'odir', help='Output directory', default='./')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    fastq_screen_list = build_filelist(args.fastq_screen_files, args.fastq_screen_list)
    fastq_screen_report = create_fastq_screen_report(fastq_screen_list)
    save_fastq_screen_report(fastq_screen_report, args.project_name, args.odir)


if __name__=="__main__":
    main()

    
