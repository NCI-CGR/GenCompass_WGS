#!/usr/bin/env python3

import json
import pandas as pd 
import argparse
import os
import datetime



def create_cmm_report(json_filename):
    with open(json_filename, "r") as fastp_json:
        json_data = json.load(fastp_json)
    alignment = json_data['report_saved_raw_data']['multiqc_picard_AlignmentSummaryMetrics']
    alignment = pd.DataFrame(alignment).T.reset_index().rename(columns={'index': 'SampleName'})

    insert_size = json_data['report_saved_raw_data']['multiqc_picard_insertSize']
    insert_size = pd.DataFrame(insert_size).T.reset_index()
    insert_size = insert_size[['SAMPLE_NAME', 'MEDIAN_INSERT_SIZE']]
    insert_size.rename(columns = {'SAMPLE_NAME': 'SampleName'}, inplace=True)


    bias = json_data['report_saved_raw_data']['multiqc_picard_gcbias']
    bias = pd.DataFrame(bias).T.reset_index().rename(columns={'index': 'SampleName'})
    bias = bias[['SampleName', 'AT_DROPOUT', 'GC_DROPOUT']]

    quality_yield = json_data['report_saved_raw_data']['multiqc_picard_QualityYieldMetrics']
    quality_yield = pd.DataFrame(quality_yield).T.reset_index().rename(columns={'index': 'SampleName'})
    quality_yield = quality_yield[['SampleName', 'TOTAL_BASES', 'Q20_BASES', 'Q30_BASES']]  

    cmm_report = (
        pd
            .merge(alignment, insert_size, on='SampleName', how = 'outer')
            .merge(bias, on='SampleName', how = 'outer')
            .merge(quality_yield, on='SampleName', how = 'outer')
    )
    return cmm_report


def save_report(cmm_report, bammetrics_report, odir, project):
    current_date = datetime.datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"
    cmm_report_filename = os.path.join(odir, f'{project}.post_mapping_report-{date_suffix}.xlsx')
    # shutil.copy(bammetrics_report, cmm_report_filename)

    thousands_columns = ['TOTAL_READS', 'PF_READS',
       'PF_NOISE_READS', 'PF_READS_ALIGNED', 
       'PF_ALIGNED_BASES', 'PF_HQ_ALIGNED_READS', 'PF_HQ_ALIGNED_BASES',
       'PF_HQ_ALIGNED_Q20_BASES', 'PF_HQ_MEDIAN_MISMATCHES',
       'PF_MISMATCH_RATE', 'PF_HQ_ERROR_RATE', 'PF_INDEL_RATE',
       'MEAN_READ_LENGTH', 'READS_ALIGNED_IN_PAIRS',
       'PF_READS_IMPROPER_PAIRS',
       'MEDIAN_INSERT_SIZE',
       'TOTAL_BASES', 'Q20_BASES', 'Q30_BASES']
    percent_columns = ['PCT_PF_READS','PCT_PF_READS_ALIGNED',
       'PCT_READS_ALIGNED_IN_PAIRS', 
       'PCT_PF_READS_IMPROPER_PAIRS',
       'PCT_CHIMERAS', 'PCT_ADAPTER']

    bammetrics_report_data = pd.read_csv(bammetrics_report, sep='\t')

    bammetrics_column_descriptions = [
        ['Sample', 'Sample name'],
        ['GENOME_TERRITORY', 'The number of non-N bases in the genome reference over which coverage will be evaluated.'],
        ['MEAN_COVERAGE', 'The mean coverage in bases of the genome territory, after all filters are applied.'],
        ['SD_COVERAGE', 'The standard deviation of coverage of the genome after all filters are applied.'],
        ['MEDIAN_COVERAGE', 'The median coverage in bases of the genome territory, after all filters are applied.'], 
        ['MAD_COVERAGE', 'The median absolute deviation of coverage of the genome after all filters are applied.'],
        ['PCT_EXC_MAPQ', 'The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20).'],
        ['PCT_EXC_DUPE','The fraction of aligned bases that were filtered out because they were in reads marked as duplicates.'],
        ['PCT_EXC_UNPAIRED', 'The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair.'],
        ['PCT_EXC_BASEQ','The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20).'],
        ['PCT_EXC_OVERLAP', 'The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.'],
        ['PCT_EXC_CAPPED', 'The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x).'],
        ['PCT_EXC_TOTAL', 'The total fraction of aligned bases excluded due to all filters.'], 
        ['PCT_1X', 'The fraction of bases that attained at least 1X sequence coverage in post-filtering bases.'],
        ['PCT_5X', 'The fraction of bases that attained at least 5X sequence coverage in post-filtering bases.'], 
        ['PCT_10X', 'The fraction of bases that attained at least 10X sequence coverage in post-filtering bases.'],
        ['PCT_15X', 'The fraction of bases that attained at least 15X sequence coverage in post-filtering bases.'],
        ['PCT_20X', 'The fraction of bases that attained at least 20X sequence coverage in post-filtering bases.'],
        ['PCT_25X', 'The fraction of bases that attained at least 25X sequence coverage in post-filtering bases.'],
        ['PCT_30X', 'The fraction of bases that attained at least 30X sequence coverage in post-filtering bases.'],
        ['PCT_40X', 'The fraction of bases that attained at least 40X sequence coverage in post-filtering bases.'],
        ['PCT_50X', 'The fraction of bases that attained at least 50X sequence coverage in post-filtering bases.'],
        ['PCT_60X', 'The fraction of bases that attained at least 60X sequence coverage in post-filtering bases.'],
        ['PCT_70X', 'The fraction of bases that attained at least 70X sequence coverage in post-filtering bases.'],
        ['PCT_80X', 'The fraction of bases that attained at least 80X sequence coverage in post-filtering bases.'],
        ['PCT_90X', 'The fraction of bases that attained at least 90X sequence coverage in post-filtering bases.'], 
        ['PCT_100X', 'The fraction of bases that attained at least 100X sequence coverage in post-filtering bases.'],
        ['HET_SNP_SENSITIVITY', 'The theoretical HET SNP sensitivity.'],
        ['HET_SNP_Q' 'The Phred Scaled Q Score of the theoretical HET SNP sensitivity.'],
        

    ]

    bammetrics_column_descriptions = pd.DataFrame(bammetrics_column_descriptions, columns=['Column', 'Description'])

    
    with pd.ExcelWriter(cmm_report_filename, engine_kwargs={'options':{'strings_to_formulas': False}}) as writer:
        bammetrics_report_data.to_excel(writer, sheet_name='bammetrics report', index=False)
        cmm_report.to_excel(writer, sheet_name="cmm report", index=False)
        bammetrics_column_descriptions.to_excel(writer, sheet_name='bammetrics column descriptions', index=False)
        
        workbook  = writer.book
        # thousands_format = workbook.add_format({'num_format': '#,##0'})
        thousands_format = workbook.add_format({'num_format': '#,###'})
        percent_format = workbook.add_format({'num_format': '0.00%'})
        
        
        # format bammetrics report
        worksheet = workbook.get_worksheet_by_name('bammetrics report')
        index_no = bammetrics_report_data.columns.get_loc('GENOME_TERRITORY')
        worksheet.set_column('B:B', None, thousands_format)
        for col in bammetrics_report_data.columns:
            if col.startswith('PCT_'):
                index_no = bammetrics_report_data.columns.get_loc(col)
                worksheet.set_column(index_no, index_no, None, percent_format)



        # format collectmultiplemetrics report
        worksheet = workbook.get_worksheet_by_name('cmm report')
        thousands_format = workbook.add_format({'num_format': '#,##0'})
        
        for col in thousands_columns:
            index_no = cmm_report.columns.get_loc(col)
            worksheet.set_column(index_no, index_no, None, thousands_format)
        
        
        for col in percent_columns:
            index_no = cmm_report.columns.get_loc(col)
            worksheet.set_column(index_no, index_no, None, percent_format)


        
        #     elif cell.value in percent_columns:
        #         cell.style = 'Percent'




def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cmm_multiqc_json', help='MultiQC json data generated using Parabricks CollectMultipleMetrics data')
    parser.add_argument('--bammetrics_report', help='bammetrics report')
    parser.add_argument('-odir', default='./')
    parser.add_argument('--project', help='project name')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    cmm_report = create_cmm_report(args.cmm_multiqc_json)
    
    save_report(cmm_report, args.bammetrics_report, args.odir, args.project)


if __name__=="__main__":
    main()
    
    
