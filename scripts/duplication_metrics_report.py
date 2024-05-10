#!/usr/bin/env python3

import pandas as pd
import click
import os

def read_fq2bam_metrics_file(fq2bam_metrics_file):
    skiprows = 0
    for line in open(fq2bam_metrics_file):
        skiprows = skiprows + 1
        if line.startswith('## METRICS CLASS'):
            break
    histograme_line_number = 0
    for line in open(fq2bam_metrics_file):
        if line.startswith('## HISTOGRAM'):
            break
        histograme_line_number = histograme_line_number + 1
    nrows = histograme_line_number - skiprows - 2
    sample_id = os.path.basename(fq2bam_metrics_file).split('.')[0]
    sample_metrics = pd.read_csv(fq2bam_metrics_file, skiprows=skiprows, nrows=nrows, sep = '\t')
    sample_metrics.insert(0, 'SAMPLE_ID', sample_id)
    return sample_metrics


def calculate_sample_level_metrics(library_metrics):
    def calculate_percent_duplication(row):
        percent_duplication = (row['UNPAIRED_READ_DUPLICATES'] + row['READ_PAIR_DUPLICATES'] * 2) / (row['UNPAIRED_READS_EXAMINED'] + row['READ_PAIRS_EXAMINED'] * 2)
        return percent_duplication
    sample_metrics = library_metrics.groupby('SAMPLE_ID').sum(numeric_only=True).reset_index()
    sample_metrics['PERCENT_DUPLICATION'] = sample_metrics.apply(calculate_percent_duplication, axis=1)
    return sample_metrics


@click.command()
@click.option('--fq2bam_metrics_fof')
@click.option('--name')
@click.option('-odir', '--output_directory')
def create_duplication_report_tabel(fq2bam_metrics_fof, name, output_directory):
    metrics_report = [read_fq2bam_metrics_file(file.strip()) for file in open(fq2bam_metrics_fof)]
    metrics_report = pd.concat(metrics_report)

    # metrics_report = metrics_report[['SAMPLE_ID', 'LIBRARY', 'READ_PAIR_DUPLICATES', 'READ_PAIR_OPTICAL_DUPLICATES', 'PERCENT_DUPLICATION']]
    
    column_descriptions = pd.DataFrame(
        [
            ["SAMPLE_ID", "Unique CGR Sample ID"],
            ["LIBRARY", "Library"],
            ["READ_PAIR_DUPLICATES", "The number of read pairs that were marked as duplicates."],
            ["READ_PAIR_OPTICAL_DUPLICATES", "The number of read pairs duplicates that were caused by optical duplication. Value is always < READ_PAIR_DUPLICATES, which counts all duplicates regardless of source."],
            ["PERCENT_DUPLICATION" , "The fraction of mapped sequence that is marked as duplicate."]
        ], 
        columns = ['Metric', 'Definition']
        
    )


    output_filename = os.path.join(output_directory, f'{name}_mapping_duplication_metrics_library.xlsx')
    with pd.ExcelWriter(output_filename) as writer:
        metrics_report.to_excel(writer, sheet_name='Duplication Metrics Report', index=False)
        column_descriptions.to_excel(writer, sheet_name="Column Descriptions", index=False)
    
    sample_metrics_report = calculate_sample_level_metrics(metrics_report)
    output_filename = os.path.join(output_directory, f'{name}_mapping_duplication_metrics_sample.xlsx')
    with pd.ExcelWriter(output_filename) as writer:
        sample_metrics_report.to_excel(writer, sheet_name='Duplication Metrics Report', index=False)
        column_descriptions.to_excel(writer, sheet_name="Column Descriptions", index=False)

if __name__ == "__main__":
    create_duplication_report_tabel()