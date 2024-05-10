#!/usr/bin/env python3

import pandas as pd
import numpy as np


def read_bedfile(bedfile: str) -> pd.DataFrame:
    """Reads interval bedfile into a pandas dataframe

    Args:
        bedfile (str): path to bedfile

    Returns:
        pd.DataFrame: bed dataframe
    """
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
        for bnum in range(bin_number + 1):
            if bin_totals.loc[bnum] < nucleotides_per_bin * 1/2:
                lower_bin_size = bin_totals.loc[bnum - 1] if bnum > 0 else np.inf
                upper_bin_size = bin_totals.loc[bnum + 1] if bnum < bin_number else np.inf
                bin_delta = -1 if lower_bin_size < upper_bin_size else 1
                for index, row in chr_intervals[chr_intervals['bin'] == bnum].iterrows():
                    chr_intervals.loc[index,'bin'] = bnum + bin_delta
    return chr_intervals

def get_interval_names(bedfile, nucleotides_per_bin):
    intervals = read_bedfile(bedfile)
    chromosomes = intervals['chromosome'].unique()
    interval_names = []
    for chr in chromosomes:
        chromosome_interval = bin_intervals(intervals=intervals, chromosome=chr, nucleotides_per_bin=nucleotides_per_bin)
        bins = chromosome_interval['bin'].unique()
        num_bins = len(bins)
        files_created = []
        for bin_id, df_bin_id in zip(range(num_bins), bins):
            if num_bins > 1:
                interval_names.append(f'{chr}.{bin_id}')
            else:
                interval_names.append(f'{chr}')
    return interval_names
