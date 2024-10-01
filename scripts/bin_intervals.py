#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import argparse
import logging
logging.basicConfig(encoding='utf-8', level=logging.DEBUG)


def read_bedfile(bedfile: str) -> pd.DataFrame:
    """Reads interval bedfile into a pandas dataframe

    Args:
        bedfile (str): path to bedfile

    Returns:
        pd.DataFrame: bed dataframe
    """
    logging.info('Reading bed file')
    with open(bedfile) as inbed:
        interval_list = []
        for line in inbed:
            splitline = line.strip('\n').split('\t')
            if len(splitline) >= 3 and not line.startswith('@'):
                interval_list.append(splitline[:3])

    intervals = pd.DataFrame(interval_list, columns=['chromosome', 'start', 'stop'])
    intervals['start'] = pd.to_numeric(intervals['start'], downcast='integer')
    intervals['stop'] = pd.to_numeric(intervals['stop'], downcast='integer')
    return intervals

def bin_intervals(intervals:pd.DataFrame, chromosome:str, nucleotides_per_bin:int=85000000) -> pd.DataFrame:
    """Bins the bedfile such that a bin attempts to hold close to a set number of nucleotides per interval bin.


    Args:
        intervals (pd.DataFrame): interval dataframe to bin
        chromosome (str): chromosome to perform binning on
        nucleotides_per_bin (int, optional): Number of nucleotides each bin should attempt to hold. Defaults to 85000000.

    Returns:
        pd.DataFrame: interval file with bins included as a column
    """
    logging.info(f'Binning intervals : {chromosome}')
    chr_intervals = intervals[intervals['chromosome'] == chromosome].copy()
    chr_intervals['size'] = chr_intervals['stop'] - chr_intervals['start']

    bin_number = 0
    bin_size = 0
    
    for index, row in chr_intervals.iterrows():
        
        if bin_size > nucleotides_per_bin or ((bin_size + row['size']) > nucleotides_per_bin * 3/4 and bin_size > 0):
            bin_number = bin_number + 1
            bin_size = 0
        chr_intervals.loc[index,'bin'] = bin_number
        bin_size = bin_size + row['size']
    
    bin_totals = chr_intervals.groupby('bin')['size'].sum()
    
    if bin_number > 0:
        logging.info('Redistributing small bins')
        for bnum in range(bin_number + 1):
            if bin_totals.loc[bnum] < nucleotides_per_bin * 1/2:
                lower_bin_size = bin_totals.loc[bnum - 1] if bnum > 0 else np.inf
                upper_bin_size = bin_totals.loc[bnum + 1] if bnum < bin_number else np.inf
                bin_delta = -1 if lower_bin_size < upper_bin_size else 1
                for index, row in chr_intervals[chr_intervals['bin'] == bnum].iterrows():
                    chr_intervals.loc[index,'bin'] = bnum + bin_delta
    return chr_intervals

def save_interval_bins(chromosome_interval:pd.DataFrame, odir:str, chromosome:str)->list:
    """Saves the intervals into individual bed files based on what bin was generated.

    Args:
        chromosome_interval (pd.DataFrame): interval file with bins as a column
        odir (str): output directory
        chromosome (str): chromosome of the interval file

    Returns:
        list: Files created
    """
    logging.info(f'Saving binned bed files: {chromosome}')
    bins = chromosome_interval['bin'].unique()
    num_bins = len(bins)
    files_created = []
    for bin_id, df_bin_id in zip(range(num_bins), bins):
        if num_bins > 1:
            ofile = os.path.join(odir, f'{chromosome}.{bin_id}.bed')
        else:
            ofile = os.path.join(odir, f'{chromosome}.bed')
        chr_bin_data: pd.DataFrame = chromosome_interval[chromosome_interval['bin'] == df_bin_id]
        chr_bin_data = chr_bin_data[['chromosome', 'start', 'stop']]
        chr_bin_data.to_csv(ofile, sep='\t', index=False, header=False)
        files_created.append(ofile)
    return files_created


    
def parse_args():
    parser= argparse.ArgumentParser(description='Bin bed file into separate bed files, where each bin attempts to hold a set number of nucleotides')
    parser.add_argument('-b','--bedfile', action='store', dest='bedfile')
    parser.add_argument('-n','--nucleotides_per_bin', action='store', dest='nucleotides_per_bin',default=85000000, type=int,  help='desired number of nucleotides per bin.' )
    parser.add_argument('-odir', '--output_directory', action='store', dest='odir')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    logging.basicConfig(filename=os.path.join(args.odir, 'bin_intervals.log'), encoding='utf-8', level=logging.DEBUG)
    os.makedirs(args.odir, exist_ok=True)
    intervals = read_bedfile(args.bedfile)
    chromosomes = intervals['chromosome'].unique()
    files_created = []
    for chr in chromosomes:
        chr_split = bin_intervals(intervals=intervals, chromosome=chr, nucleotides_per_bin=args.nucleotides_per_bin)
        chr_files_created = save_interval_bins(chr_split, args.odir, chr)
        files_created = files_created + chr_files_created
    with open(os.path.join(args.odir, 'binned_intervals_FOF.txt'), 'w') as ofile:
        ofile.write('\n'.join(files_created))

if __name__=="__main__":
    main()
