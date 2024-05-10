#!/usr/bin/env python3

import argparse 
import json
import os
import pandas as pd
import datetime
from helpers import build_filelist, Manifest

    
def parse_fastp_json(json_filename):
    """Parse single json file into dict of column name:value pairs

    Args:
        json_filename (str): fastp json filename

    Returns:
        dict: parsed fastp data
    """
    sample_name = {"sample": os.path.basename(json_filename).strip('.json')}
    
    with open(json_filename, "r") as fastp_json:
        json_data = json.load(fastp_json)

        ## Define groups of metrics to include in output (dictionaries)
        before_filtering = json_data['summary']['before_filtering']
        after_filtering = json_data['summary']['after_filtering']
        filtering_result = json_data['filtering_result']
        dup_rate = json_data['duplication']
        insert_size = json_data['insert_size']
        adapter_cutting = json_data['adapter_cutting']

        #count total adapters cut 
        adapter_cutting_results = {}
        adapter_cutting_results['Num_Adapters_Removed_R1'] = sum(adapter_cutting['read1_adapter_counts'].values())
        adapter_cutting_results['Num_Adapters_Removed_R2'] = sum(adapter_cutting['read2_adapter_counts'].values())

        
        ## Delete insert size histogram
        del insert_size['histogram']
        
        ## Edit some of the dictionary keys to include more information about the metrics
        def add_key_prefix(dictionary, prefix):
            new_dictionary = {prefix+k: v for k, v in dictionary.items()}
            return new_dictionary
        
        before_filtering = add_key_prefix(before_filtering, "before_filtering_")
        after_filtering = add_key_prefix(after_filtering, "after_filtering_")
        dup_rate = add_key_prefix(dup_rate, "duplication_")
        insert_size = add_key_prefix(insert_size, "insert_size_")
        
        ## Merge dictionaries to consolidate data for tsv
        tsv_data = {**sample_name, **before_filtering, **after_filtering, **filtering_result, **dup_rate, **insert_size, **adapter_cutting_results}
    
    return tsv_data

# def read_keytable(keytable_filename):
#     keytable = pd.read_csv(keytable_filename)
#     keytable.rename(columns={'SampleID': 'USU_ID'}, inplace=True)
#     keytable['SampleID'] = keytable['Description'].apply(lambda desc: desc.split('-')[-1])
#     return keytable
    
def build_fastp_report(fastp_json_files: list, manifest: Manifest, project_name: str, sequencing_type: str):
    fastp_report = [parse_fastp_json(x) for x in fastp_json_files]
    fastp_report = pd.DataFrame(fastp_report)
    fastp_report.rename(columns = {
        'sample': 'Lane_ID',
        'Flowcell': 'FlowCell_ID',
        'before_filtering_total_reads': 'Total_Raw_Reads_Lane',
        'before_filtering_total_bases': 'R1_R2_Raw_Bases',
        'before_filtering_q30_bases': 'R1_R2_Raw_Q30_bases',
        'before_filtering_read1_mean_length': 'Untrimmed_Mean_R1_Length', 
        'before_filtering_read2_mean_length': 'Untrimmed_Mean_R2_Length', 
        'before_filtering_gc_content': 'R1_R2_Raw_GC_Content(%)',
        'after_filtering_total_reads': 'Total_Filtered_Reads_Lane',
        'after_filtering_total_bases': 'R1_R2_Filtered_Bases',
        'after_filtering_q30_bases': 'R1_R2_Filtered_Q30_bases',
        'after_filtering_read1_mean_length': 'Trimmed_Mean_R1_Length', 
        'after_filtering_read2_mean_length': 'Trimmed_Mean_R2_Length', 
        'after_filtering_gc_content': 'Trimmed_Filtered_R1_R2_GC_Content(%)',
        # 'passed_filter_reads': 'Num_Reads_Passing_QC_Filters',
        # 'low_quality_reads': 'Percent_R1_R2_less_than_Q25',
        # 'too_many_N_reads': 'Percent_R1_R2_greater_than_10_Ns',
        # 'low_complexity_reads': 'Percent_R1_R2_Low_Complexity',
        # 'too_short_reads': 'Percent_R1_R2_less_than_50_nt',
        'duplication_rate': 'Raw_R1_R2_Duplication_Rate(%)'
    }, inplace=True)
    fastp_report[manifest.sample_run_id_column] = fastp_report['Lane_ID'].apply(lambda lane: lane.split('_')[0])

    #Change to percentages
    fastp_report['Percent_R1_R2_less_than_Q25'] = fastp_report['low_quality_reads'] / fastp_report['Total_Raw_Reads_Lane'] 
    fastp_report['Percent_R1_R2_greater_than_10_Ns'] = fastp_report['too_many_N_reads'] / fastp_report['Total_Raw_Reads_Lane'] 
    fastp_report['Percent_R1_R2_Low_Complexity'] = fastp_report['low_complexity_reads'] / fastp_report['Total_Raw_Reads_Lane'] 
    fastp_report['Percent_R1_R2_less_than_50_nt'] = fastp_report['too_short_reads'] / fastp_report['Total_Raw_Reads_Lane'] 
    fastp_report['Percent_R1_R2_Passing_QC_Filters'] = fastp_report['Total_Filtered_Reads_Lane'] / fastp_report['Total_Raw_Reads_Lane'] 

    fastp_report['Untrimmed_Mean_Length (R1;R2)'] = fastp_report.apply(lambda row: f'{row["Untrimmed_Mean_R1_Length"]};{row["Untrimmed_Mean_R2_Length"]}' ,axis=1)
    fastp_report['Trimmed_Mean_Length (R1;R2)'] = fastp_report.apply(lambda row: f'{row["Trimmed_Mean_R1_Length"]};{row["Trimmed_Mean_R2_Length"]}' ,axis=1)

    fastp_report = fastp_report.merge(manifest.manifest, on = manifest.sample_run_id_column, how = 'left')
    fastp_report.rename(columns = {
        'Flowcell': 'FlowCell_ID',
        manifest.sample_run_id_column: 'Sample_Run_ID',
        manifest.sample_id_column: 'Sample_ID'
    }, inplace=True)
    fastp_report['Project_Name'] = project_name
    fastp_report['Sequencing_Type'] = sequencing_type

    fastp_report['Total_Raw_Read_Pairs_Lane'] = fastp_report['Total_Raw_Reads_Lane']  // 2
    fastp_report['Total_Filtered_Read_Pairs_Lane'] = fastp_report['Total_Filtered_Reads_Lane']  // 2
    sample_summary = (
        fastp_report[['Sample_ID', 'Total_Raw_Reads_Lane', 'Total_Filtered_Reads_Lane', 'Total_Raw_Read_Pairs_Lane', 'Total_Filtered_Read_Pairs_Lane']].groupby('Sample_ID')
            .sum().reset_index()
            .rename(columns = {
                'Total_Raw_Reads_Lane': 'Total_Raw_Reads_Sample', 
                'Total_Filtered_Reads_Lane': 'Total_Filtered_Reads_Sample',
                'Total_Raw_Read_Pairs_Lane': 'Total_Raw_Read_Pairs_Sample',
                'Total_Filtered_Read_Pairs_Lane': 'Total_Filtered_Read_Pairs_Sample'

                })
    )

    fastp_report = fastp_report.merge(sample_summary, on = 'Sample_ID')

    fastp_report = fastp_report[[ 
        'Sample_ID',
        'FlowCell_ID',
        'Sample_Run_ID',
        'Lane_ID',
        'Project_Name',
        'Sequencing_Type',
        # 'Total_Raw_Reads_Lane',
        # 'Total_Raw_Reads_Sample',
        'Total_Raw_Read_Pairs_Lane',
        'Total_Raw_Read_Pairs_Sample',
        'R1_R2_Raw_Bases',
        'R1_R2_Raw_Q30_bases',
        'Untrimmed_Mean_Length (R1;R2)',
        'R1_R2_Raw_GC_Content(%)',
        'Raw_R1_R2_Duplication_Rate(%)',
        # 'Total_Filtered_Reads_Lane',
        # 'Total_Filtered_Reads_Sample',
        'Total_Filtered_Read_Pairs_Lane',
        'Total_Filtered_Read_Pairs_Sample',
        'R1_R2_Filtered_Bases',
        'R1_R2_Filtered_Q30_bases',
        'Trimmed_Mean_Length (R1;R2)',
        'Trimmed_Filtered_R1_R2_GC_Content(%)',
        'Percent_R1_R2_Passing_QC_Filters',
        'Percent_R1_R2_less_than_Q25',
        'Percent_R1_R2_greater_than_10_Ns',
        'Percent_R1_R2_Low_Complexity',
        'Percent_R1_R2_less_than_50_nt',
        'Num_Adapters_Removed_R1',
        'Num_Adapters_Removed_R2'
    ]]

    
    return fastp_report

    
    
def save_fastp_report(fastp_report, odir, project_name):
    # Create and save column descriptions
    column_descriptions = pd.DataFrame(
        [
            ['Sample_ID', 'Unique CGR Sample ID'],
            ['FlowCell_ID', 'Flowcell(s) for a Sample'],
            ['Sample_Run_ID', 'Unique USU ID, Generally one USU ID corresponds to one CGR Sample ID unless there are top-offs'],
            ['Lane_ID', 'The Flowcell Lane IDs for a sample (for most samples, there are 4 lanes/sample)'],
            ['Project_Name', 'Project Name'],
            ['Sequencing_Type', 'Sequencing Type'],
            ['Total_Raw_Read_Pairs_Lane', 'Number of total raw read pairs for a lane'], 
            ['Total_Raw_Read_Pairs_Sample', 'Number of total raw read pairs for a sample'],
            ['R1_R2_Raw_Bases', 'Number of bases in the raw R1 and R2 reads combined together'],
            ['R1_R2_Raw_Q30_bases', 'Number of Q30 bases in the raw R1 and R2 reads combined together'],
            ['Untrimmed_Mean_Length (R1;R2)', 'The mean read length of raw untrimmed R1 and R2 reads'],
            ['R1_R2_Raw_GC_Content(%)', 'The GC content of raw R1 and R2 reads expressed as percentage'],
            ['Raw_R1_R2_Duplication_Rate(%)', 'The duplication rate of raw R1 and R2 reads expressed as percentage'],
            ['Total_Filtered_Read_Pairs_Lane', 'Filtered read pairs/lane that pass the criteria:  base quality > 25 in at least 40% of the bases/read pair, number of N\'s < 10/read pair, average quality > 25/read pair, minimum read length = 50, reads pairs > 30% complexity'],
            ['Total_Filtered_Read_Pairs_Sample', 'Filtered read pairs/sample that pass the criteria as above'],
            ['R1_R2_Filtered_Bases', 'Number of bases in the filtered R1 and R2 reads combined together'],
            ['R1_R2_Filtered_Q30_bases', 'Number of Q30 bases in the filtered R1 and R2 reads combined together'],
            ['Trimmed_Mean_Length (R1;R2)', 'The mean read length of raw trimmed R1 and R2 reads'],
            ['Trimmed_Filtered_R1_R2_GC_Content(%)', 'The GC content of trimmed & filtered R1 and R2 reads expressed as percentage'],
            ['Percent_R1_R2_Passing_QC_Filters', '% paired-end reads passing the QC filters'],
            ['Percent_R1_R2_less_than_Q25', '% paired-end reads less than Q25'],
            ['Percent_R1_R2_greater_than_10_Ns', '% paired-end reads with 10 or more Ns'],
            ['Percent_R1_R2_Low_Complexity', '% paired-end reads with less than 30% complexity'],
            ['Percent_R1_R2_less_than_50_nt', '% paired-end reads that are less than 50 bases'],
            ['Num_Adapters_Removed_R1', 'Number of adapters trimmed from R1 reads'],
            ['Num_Adapters_Removed_R2', 'Number of adapters trimmed from R2 reads']
        ],
        columns=['Header', 'Description']
    )
    
    current_date = datetime.datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"
    with pd.ExcelWriter(os.path.join(odir, f'{project_name}.pre-mapping-QC-A-fastp_report-{date_suffix}.xlsx')) as writer:  
        fastp_report.to_excel(writer, sheet_name='Fastp Report', index=False)
        column_descriptions.to_excel(writer, sheet_name='column-header-descriptions', index=False)
        workbook  = writer.book
        # worksheet = writer.sheets['Fastp_Report']
        worksheet = workbook.get_worksheet_by_name('Fastp Report')
        # format read as with comma between thousands
        thousands_format = workbook.add_format({'num_format': '#,##0'})
        thousands_unit_columns = [
            # 'Total_Raw_Reads_Lane',
            # 'Total_Raw_Reads_Sample',
            'Total_Raw_Read_Pairs_Lane',
            'Total_Raw_Read_Pairs_Sample',
            'R1_R2_Raw_Bases',
            'R1_R2_Raw_Q30_bases',
            # 'Total_Filtered_Reads_Lane',
            # 'Total_Filtered_Reads_Sample',
            'Total_Filtered_Read_Pairs_Lane',
            'Total_Filtered_Read_Pairs_Sample',
            'R1_R2_Filtered_Bases',
            'R1_R2_Filtered_Q30_bases',
            'Num_Adapters_Removed_R1',
            'Num_Adapters_Removed_R2']
        for col in thousands_unit_columns:
                index_no = fastp_report.columns.get_loc(col)
                worksheet.set_column(index_no, index_no, None, thousands_format)

        # format percentages
        percent_format = workbook.add_format({'num_format': '0.00%'})
        percent_unit_columns = [
            'Percent_R1_R2_Passing_QC_Filters',
            'Percent_R1_R2_less_than_Q25',
            'Percent_R1_R2_greater_than_10_Ns',
            'Percent_R1_R2_Low_Complexity',
            'Percent_R1_R2_less_than_50_nt',
            'R1_R2_Raw_GC_Content(%)', 
            'Raw_R1_R2_Duplication_Rate(%)', 
            'Trimmed_Filtered_R1_R2_GC_Content(%)']
        for col in percent_unit_columns:
                index_no = fastp_report.columns.get_loc(col)
                worksheet.set_column(index_no, index_no, None, percent_format)

        

    
    


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fastp_json', nargs='*', required=False,  help='Fastp JSON files to combine into single file', )
    parser.add_argument('-l', '--fastp_json_list', required=False, help='List of fastp JSON files to combine into single file')
    parser.add_argument('-m', '--manifest', help='Manifest file containing Sample ID')
    parser.add_argument('-p', '--project_name', help='Name of project')
    parser.add_argument('-s', '--sequencing_type', help='Sequencing done on samples')
    parser.add_argument('-odir','--output_directory', dest='odir', default = './', help='Output directory')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    fastp_json_list = build_filelist(args.fastp_json, args.fastp_json_list)

    manifest = Manifest(args.manifest)
    fastp_report = build_fastp_report(fastp_json_list, manifest, args.project_name, args.sequencing_type)
    
    save_fastp_report(fastp_report, args.odir, args.project_name)

if __name__=="__main__":
    main()
