#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import datetime

def rename_sample(samplename):
    new_name = samplename[:len(samplename)//2]
    return new_name

def infer_sex(row):
    if row['Y_depth_mean'] == 0 and row['X_het'] > 0:
        return 'Female'
    elif row['Y_depth_mean'] > 0 and row['X_het'] ==0:
        return 'Male'
    else:
        return 'Unknown'


def create_sex_report(somalier_sample_report) -> pd.DataFrame:
    somalier_sample_report.drop(columns=['paternal_id', 'maternal_id', 'sex', 'phenotype', 'original_pedigree_sex'], inplace=True)
    somalier_sample_report['Inferred_Sex'] = somalier_sample_report.apply(infer_sex, axis=1)
    return somalier_sample_report

def create_relatedness_report(somalier_relatedness_report, ibs2_cutoff=10000, relatedness_cutoff=0.1) -> pd.DataFrame:
    somalier_relatedness_report = somalier_relatedness_report[
        (somalier_relatedness_report['relatedness'] > relatedness_cutoff) & 
        (somalier_relatedness_report['ibs2'] > ibs2_cutoff)]
    return somalier_relatedness_report

def parse_args():
    parser = argparse.ArgumentParser() 
    parser.add_argument('--somalier_samples', help='somalier samples report')
    parser.add_argument('--somalier_pairs', help='somalier pairs report')
    parser.add_argument('-odir', default='./')
    parser.add_argument('--project', help='project name')
    parser.add_argument('--ibs2_cutoff', type=float, default=10000, help='IBS2 cutoff used for filtering relatedness report')
    parser.add_argument('--relatedness_cutoff', type=float, default = 0.1, help='relatedness cutoff used for filtering relatedness report')
    args = parser.parse_args()
    return args

def save_relatedness_report(relatedness_report:pd.DataFrame, odir: str, project: str):
    column_descriptions = pd.DataFrame(
        [
            ["sample_a", "First sample in pair"],
            ["sample_b", "Second sample in pair"],
            ["relatedness",	"Relatedness value between 0 and 1. By default, only sample pairs with relatedness >0.1 are printed."],
            ["ibs0", "The number of sites where one sample is hom-ref and another is hom-alt and the 2 samples shared no alleles (should approach 0 for parent-child pairs)"],
            ["ibs2", "The number of sites at which the 2 samples are both hom-ref, both het, or both hom-alt. By default, only sample pairs with IBS2>10000 are printed"],
            ["hom_concordance", "Homozygous concordance"],
            ["hets_a", "The number of sites at which sample_a was heterozygous"],
            ["hets_b", "The number of sites at which sample_b was heterozygous"],
            ["hets_ab", "Count of sites that are het in either sample and not unknown in the other sample (https://github.com/brentp/somalier/blob/master/src/somalierpkg/bitset.nim)"],
            ["shared_hets", "The number of sites at which both samples were hets"],
            ["hom_alts_a", "The number of sites at which sample_a was homozygous"],
            ["hom_alts_b", "The number of sites at which sample_b was homozygous"],
            ["shared_hom_alts", "the number of sites where both samples are homozygous alternate"],
            ["n", "Number of 0/0, 0/1 and 1/1 sites"],
            ["x_ibs0", "Number of 0/1 sites on X chromosome"],
            ["x_ibs2", "Number of 1/1 sites on X chromosome"],
            ["expected_relatedness", "Ignore, always -1"],
        ],
        columns=['Column', 'Description']
    )
    current_date = datetime.datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"
    relate_filename = os.path.join(odir, f'{project}.bam_QC_somalier_relatedness-{date_suffix}.xlsx')
    with pd.ExcelWriter(relate_filename) as writer:
        sheet_name = "somalier_relatedness"
        relatedness_report.to_excel(writer, sheet_name=sheet_name, index=False)
        column_descriptions.to_excel(writer, sheet_name="Column Descriptions", index=False)

def save_sex_report(sex_report: pd.DataFrame, odir: str, project: str):
    column_descriptions = pd.DataFrame(
        [
            ["family_id", "Family ID (same as sample ID if no family structure is provided)"],
            ["sample_id", "Sample ID"],
            ["gt_depth_mean", "Mean depth of genotyped sites"],
            ["gt_depth_sd", "Standard deviation of depth of genotyped sites"],
            ["depth_mean", "Mean depth of all sites"],
            ["depth_sd", "Standard deviation of depth of all sites"],
            ["ab_mean", "Mean allele balance"],
            ["ab_std", "Standard deviation of allele balance"],
            ["n_hom_ref", "Number of 0/0 sites"],
            ["n_het", "Number of 0/1 sites"],
            ["n_hom_alt", "Number of 1/1 sites"],
            ["n_unknown", "Number of unknown sites"],
            ["p_middling_ab", "Proportion sites with AB<0.1 or AB>0.9"],
            ["X_depth_mean", "Scaled mean depth on chrX"],
            ["X_n", "Total number of sites on chrX"],
            ["X_hom_ref", "Number of 0/0 sites on chrX"],
            ["X_het", "Number of 0/1 sites on chrX"],
            ["X_hom_alt", "Number of 1/1 sites on chrX"],
            ["Y_depth_mean", "Scaled mean depth on chrY"],
            ["Y_n", "Total number of sites on chrY"],
            ["Inferred_Sex", "Female if Y_depth_mean=0 and X_het>0. Male if Y_depth_mean>0 and X_het=0."],
        ]
    )
    current_date = datetime.datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"
    sex_filename = os.path.join(odir, f'{project}.bam_QC_somalier_sex_report-{date_suffix}.xlsx')
    with pd.ExcelWriter(sex_filename) as writer:
        sheet_name = "somalier_sex_report"
        sex_report.to_excel(writer, sheet_name=sheet_name, index=False)
        column_descriptions.to_excel(writer, sheet_name="Column Descriptions", index=False)

def main():
    args = parse_args()
    somalier_sample_report = pd.read_csv(args.somalier_samples, sep='\t')
    somalier_pairs_report = pd.read_csv(args.somalier_pairs, sep='\t')

    somalier_sex_report = create_sex_report(somalier_sample_report)
    somalier_relatedness_report = create_relatedness_report(somalier_pairs_report, ibs2_cutoff=args.ibs2_cutoff, relatedness_cutoff=args.relatedness_cutoff)

    sex_filename = os.path.join(args.odir, f'{args.project}.somalier_sex_report.tsv')
    somalier_sex_report.to_csv(sex_filename, sep='\t', index=False)
    
    save_relatedness_report(somalier_relatedness_report, args.odir, args.project)
    save_sex_report(somalier_sex_report, args.odir, args.project)


if __name__=="__main__":
    main()
