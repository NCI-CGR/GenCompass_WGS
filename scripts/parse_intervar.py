#!/usr/bin/env python3

import argparse 
import pandas as pd
import numpy as np
import re 

def extract_significance(intervar_info):
    significance = intervar_info.strip('InterVar: ').split('PVS1')[0].strip()
    return significance

def extract_PVS1(intervar_info):
    match = re.search('PVS1=(\d+)', intervar_info)
    if match:
        pvs1 = match.group(1)
        return float(pvs1)
    else:
        return np.nan

def parse_intervar(intervar_filename):
    intervar = pd.read_csv(intervar_filename, sep = '\t')
    sig_col = ' InterVar: InterVar and Evidence '


    intervar['InterVar.significance'] = intervar[sig_col].apply(extract_significance)
    intervar['InterVar.PVS1'] = intervar[sig_col].apply(extract_PVS1)

    intervar = intervar[['#Chr', 'Start','End', 'Ref', 'Alt', sig_col, 'InterVar.significance', 'InterVar.PVS1']]
    intervar.rename(columns = {sig_col: sig_col.strip()}, inplace=True)
    intervar.drop_duplicates(inplace=True)
    # intervar.rename(columns = {'#Chr': 'CHROM', 'Start': 'POS', 'Ref': 'REF', 'Alt': 'ALT'}, inplace=True)
    return intervar

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--intervar', action='store', dest='intervar')
    parser.add_argument('--output', action='store', dest='output')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    parsed_intervar = parse_intervar(args.intervar)
    parsed_intervar.to_csv(args.output, sep = '\t', index=False)

if __name__=="__main__":
    main()


