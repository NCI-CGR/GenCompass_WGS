#!/usr/bin/env python3
import argparse 
import os 
import pandas as pd
import datetime
from helpers import build_filelist


def read_sample_kraken_report(kraken_report_filename: str):
    sample_name = os.path.basename(kraken_report_filename).split('.')[0]
    report = pd.read_csv(kraken_report_filename, sep='\t', index_col=0, header=None, names=[sample_name])
    return report

def create_classified_reads_report(kraken_report_filenames: list) -> pd.DataFrame:
    kraken_reports = [read_sample_kraken_report(f) for f in kraken_report_filenames]
    kraken_reports = pd.concat(kraken_reports, axis=1).fillna(0)
    return kraken_reports

def build_readcount(kraken_output_filename):
    total_classified = 0
    total_unclassified = 0
    for line in open(kraken_output_filename):
        if line.startswith("C"):
            total_classified = total_classified + 1
        elif line.startswith("U"):
            total_unclassified = total_unclassified + 1
    total_processed = total_unclassified + total_classified
    return {"Total Classified": total_classified,
            "Total Unclassified": total_unclassified,
            "Total Reads Processed": total_processed}

def create_read_totals_report(output_filenames):
    read_totals = pd.DataFrame(index=["Total Reads Processed", "Total Classified", "Total Unclassified"])
    for filename in output_filenames:
        sample_name = os.path.basename(filename).split('.')[0]
        readcount_dict = build_readcount(filename)
        for key, value in readcount_dict.items():
            read_totals.at[key, sample_name] = value
        read_totals[sample_name] = read_totals[sample_name].astype(int)
    return read_totals

def create_kraken_report(kraken_report_filenames: list, kraken_output_filenames: list) -> pd.DataFrame:
    classified_reads = create_classified_reads_report(kraken_report_filenames)
    read_totals = create_read_totals_report(kraken_output_filenames)
    kraken_report = pd.concat([read_totals, classified_reads]).astype(int)
    return kraken_report

def save_kraken_report(kraken_report: pd.DataFrame, odir: str, project: str):
    row_descriptions = pd.DataFrame([
        ["SampleName", "Sample ID"],
        ["Total Reads Processed", "Total number of reads processed by Kraken (combination of all unmapped R1 and R2 reads)"],
        ["Total Classified", "Total number of reads classified by Kraken"],
        ["Total Unclassified", "Total number of reads not classified by Kraken"],
        ["Additional rows (ex. d__Viruses|f__Anelloviridae|g__Alphatorquevirus)", "Number of reads classified to the specified genus. The full taxonomy is given for each genus, with each rank denoted by a lowercase letter (ex. \"d__\" denotes the domain, \"g__\" denotes the genus)."],
    ], columns=["Row", "Description"])

    current_date = datetime.datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"

    filename = os.path.join(odir, f"{project}.bam-QC-Kraken-genus-level-{date_suffix}.xlsx")
    kraken_report = kraken_report.reset_index().rename(columns={'index': 'SampleName'})
    with pd.ExcelWriter(filename) as writer:

        kraken_report.to_excel(writer, sheet_name=f"Kraken Report", index=False)
        row_descriptions.to_excel(writer, sheet_name="Row Descriptions", index=False)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kraken_reports', nargs='+', help='Kraken MPA style reports')
    parser.add_argument('--kraken_outputs', nargs='+', help='Kraken output. Used to count number of unclassified reads')
    parser.add_argument('--project', help='name of project. used in naming output')
    parser.add_argument('-odir', '--output_directory', dest='odir', help='Output directory')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    if len(args.kraken_reports) == 1:
        kraken_report_files = build_filelist(files=None, files_fof = args.kraken_reports[0])
    else:
        kraken_report_files = args.kraken_reports
    if len(args.kraken_outputs) == 1:
        kraken_output_files = build_filelist(files=None, files_fof = args.kraken_outputs[0])
    else:
        kraken_output_files = args.kraken_outputs
    kraken_report = create_kraken_report(kraken_report_files, kraken_output_files)
    save_kraken_report(kraken_report, args.odir, args.project)

if __name__=="__main__":
    main()