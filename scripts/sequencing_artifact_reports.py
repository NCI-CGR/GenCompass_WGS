#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import datetime
from helpers import build_filelist

def read_sequencing_artifact_file(filename: str) -> pd.DataFrame:
    with open(filename) as ifile:
        skiprows = 0
        for line in ifile.readlines():
            
            skiprows = skiprows + 1
            if line.startswith('## METRICS CLASS'):
                break
    artifacts = pd.read_csv(filename, skiprows=skiprows, sep = '\t')
    return artifacts

def concat_sequencing_artifact_summaries(artifact_summaries: list) -> pd.DataFrame:
    summaries_parsed = [read_sequencing_artifact_file(filename) for filename in artifact_summaries]
    summaries_concat = pd.concat(summaries_parsed, ignore_index=True)
    return summaries_concat

def save_bait_bias_report(bait_bias_report: pd.DataFrame, odir: str, project: str):
    column_descriptions = pd.DataFrame([
        ["SAMPLE_ALIAS", "The name of the sample being assayed."],
        ["LIBRARY", "The name of the library being assayed."],
        ["REF_BASE", "The (upper-case) original base on the reference strand."],
        ["ALT_BASE", "The (upper-case) alternative base that is called as a result of DNA damage."],
        ["TOTAL_QSCORE", "The total Phred-scaled Q-score for this artifact. A lower Q-score means a higher probability that a REF_BASE:ALT_BASE observation randomly picked from the data will be due to this artifact, rather than a true variant."],
        ["WORST_CXT", "The sequence context (reference bases surrounding the locus of interest) having the lowest Q-score among all contexts for this artifact."],
        ["WORST_CXT_QSCORE", "The Q-score for the worst context."],
        ["WORST_PRE_CXT", "The pre-context (reference bases leading up to the locus of interest) with the lowest Q-score."],
        ["WORST_PRE_CXT_QSCORE", "The Q-score for the worst pre-context."],
        ["WORST_POST_CXT", "The post-context (reference bases trailing after the locus of interest) with the lowest Q-score."],
        ["WORST_POST_CXT_QSCORE", "The Q-score for the worst post-context."],
        ["ARTIFACT_NAME", "A \"nickname\" of this artifact, if it is a known error mode."],
    ])
    current_date = datetime.datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"
    ofilename = os.path.join(odir, f"{project}.bam_QC_C_sequencingArtifact.bait_bias-{date_suffix}.xlsx")
    with pd.ExcelWriter(ofilename) as writer:
        bait_bias_report.to_excel(writer, sheet_name=f"seqArt.bait", index=False)
        column_descriptions.to_excel(writer, sheet_name="Column Descriptions", index=False)



def save_pre_adapter_report(pre_adapter_report: pd.DataFrame, odir: str, project: str):
    column_descriptions = pd.DataFrame([
        ["SAMPLE_ALIAS", "The name of the sample being assayed."],
        ["LIBRARY", "The name of the library being assayed."],
        ["REF_BASE", "The (upper-case) original base on the reference strand."],
        ["ALT_BASE", "The (upper-case) alternative base that is called as a result of DNA damage."],
        ["TOTAL_QSCORE", "The total Phred-scaled Q-score for this artifact. A lower Q-score means a higher probability that a REF_BASE:ALT_BASE observation randomly picked from the data will be due to this artifact, rather than a true variant."],
        ["WORST_CXT", "The sequence context (reference bases surrounding the locus of interest) having the lowest Q-score among all contexts for this artifact."],
        ["WORST_CXT_QSCORE", "The Q-score for the worst context."],
        ["WORST_PRE_CXT", "The pre-context (reference bases leading up to the locus of interest) with the lowest Q-score."],
        ["WORST_PRE_CXT_QSCORE", "The Q-score for the worst pre-context."],
        ["WORST_POST_CXT", "The post-context (reference bases trailing after the locus of interest) with the lowest Q-score."],
        ["WORST_POST_CXT_QSCORE", "The Q-score for the worst post-context."],
        ["ARTIFACT_NAME", "A \"nickname\" of this artifact, if it is a known error mode."],
    ])
    current_date = datetime.datetime.now()
    date_suffix = f"{current_date.year}{current_date.month:02d}{current_date.day:02d}"
    ofilename = os.path.join(odir, f"{project}.bam_QC_C_sequencingArtifact.pre_adapter-{date_suffix}.xlsx")
    with pd.ExcelWriter(ofilename) as writer:
        pre_adapter_report.to_excel(writer, sheet_name=f"seqArt.pre_a", index=False)
        column_descriptions.to_excel(writer, sheet_name="Column Descriptions", index=False)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--seqart_pre', nargs='*', required=False, help='sequencingArtifact.pre_adapter_summary_metrics.txt files for samples.')
    parser.add_argument('-b', '--seqart_bait', nargs='*', required=False, help='sequencingArtifact.bait_bias_summary_metrics.txt files for samples.')
    parser.add_argument('-pl', '--seqart_pre_list', required=False, help='List of sequencingArtifact.pre_adapter_summary_metrics.txt files for samples.')
    parser.add_argument('-pb', '--seqart_bait_list', required=False, help='List of sequencingArtifact.bait_bias_summary_metrics.txt files for samples.')
    parser.add_argument('--project')
    parser.add_argument('-odir', '--output_directory', dest='odir')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    seqart_pre_filelist = build_filelist(args.seqart_pre, args.seqart_pre_list)
    seqart_bait_filelist = build_filelist(args.seqart_bait, args.seqart_bait_list)
    pre_adapter_report = concat_sequencing_artifact_summaries(seqart_pre_filelist)
    bait_bias_report = concat_sequencing_artifact_summaries(seqart_bait_filelist)

    save_pre_adapter_report(pre_adapter_report, args.odir, args.project)
    save_bait_bias_report(bait_bias_report, args.odir, args.project)
if __name__=="__main__":
    main()
