#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def classify_variant_type(ref, alt):
    """Classifiy the variant type based on length of reference and alternate nucleotides
    See https://docs.gdc.cancer.gov/Encyclopedia/pages/Variant_Type/#variant-type
    for more details on variant types

    Args:
        ref (str): reference nucleotides of variant
        alt (str): alternatve nucloetides of variant

    Returns:
        str: variant type classification
    """
    if len(alt) > len(ref):
        return 'INS'
    elif len(alt) < len(ref):
        return 'DEL'

    # Equal length
    elif len(ref) == 1 and len(alt) == 1:
        return 'SNP'
    elif len(ref) == 2 and len(alt) == 2:
        return 'DNP'
    elif len(ref) == 3 and len(alt) == 3:
        return 'TNP'
    else:
        return 'ONP'



def merge_annotations(annovar, intervar, snpeff):
    """Merge annotations into single dataframe

    Args:
        annovar (DataFrame): annovar data
        intervar (DataFrame): intervar data, relevant information extracted beforehand
        snpeff (DataFrame): SnpEff data, relevant information extracted beforehand
    """

    def get_annovar_id(row):
        # IF DELETION
        if len(row['Ref']) > len(row['Alt']):
            return f"{row['Chr']}_{row['Start']}_{row['Ref']}_{row['Alt']}"
        # IF SNP, INSERTION
        else:
            return f"{row['Chr']}_{row['End']}_{row['Ref']}_{row['Alt']}"
    def get_snpeff_id(row):
        return f"{row['CHROM']}_{row['POS']}_{row['REF']}_{row['ALT']}"    

    annovar['UID'] = annovar.apply(get_annovar_id, axis = 1)
    intervar['Chr'] = intervar['#Chr'].apply(lambda x: "chr" + str(x))
    intervar.drop(columns=['#Chr'], inplace=True)
    intervar['UID'] = intervar.apply(get_annovar_id, axis = 1)
    intervar_cols = ['UID', 'InterVar.significance', 'InterVar.PVS1']
    # intervar_cols = list(set(intervar.columns) - set(annovar.columns))
    # intervar_cols.append('UID')

    merged = pd.merge(intervar[intervar_cols], annovar, on ='UID', how = 'outer')
    merged.drop(columns = ['Chr', 'Start', 'End', 'Ref', 'Alt'], inplace=True)


    snpeff['UID'] = snpeff.apply(get_snpeff_id, axis = 1)

    merged= snpeff.merge(merged, how = 'outer', on='UID')
    merged.drop('UID', axis=1, inplace=True)
    variant_type = merged.apply(lambda row: classify_variant_type(row['REF'], row['ALT']), axis = 1)
    merged.insert(8, 'VariantType', variant_type)
    # merged['POS'] = merged['POS'].astype(str).astype(int)

    
    return merged


def save_by_feature(merged, name, odir):
    """Separate the data based on feature and save each feature set individually.

    gene: UTR, exonic, splicing, upstream, downstream
    intronic: gene intronic variants
    ncRNA: ncRNA variants
    intergenic: intergenic variants


    Args:
        merged (DataFrame): merged variant annotations
        name (str): prepended name used in file naming
        odir (str): output directory
    """
    
    merged_nona = merged[merged['Func.refGene'].notna()]
    not_ncrna = merged_nona[~merged_nona['Func.refGene'].str.contains('ncRNA')]

    protein_features = ['UTR', 'exonic','splicing', 'upstream', 'downstream']
    merged_protein = not_ncrna[not_ncrna['Func.refGene'].str.contains('|'.join(protein_features))]

    ncrna = merged_nona[merged_nona['Func.refGene'].str.contains('ncRNA')]
    intergenic = not_ncrna[not_ncrna['Func.refGene'] == 'intergenic']
    intronic = not_ncrna[not_ncrna['Func.refGene'] == 'intronic']

    merged_protein.to_csv(os.path.join(odir, f"{name}_annotated.gene.tsv"), sep = '\t', index=False)
    ncrna.to_csv(os.path.join(odir, f"{name}_annotated.ncRNA.tsv"), sep = '\t', index=False)
    intergenic.to_csv(os.path.join(odir, f"{name}_annotated.intergenic.tsv"), sep = '\t', index=False)
    intronic.to_csv(os.path.join(odir, f"{name}_annotated.intronic.tsv"), sep = '\t', index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--annovar', action='store', dest='annovar', help='Annovar output file. CSV format')
    parser.add_argument('--intervar', action='store', dest='intervar', help='InterVar output file. CSV format')
    parser.add_argument('--snpeff', action='store', dest='snpeff', help='SnpEff output, with effect extraction. TSV format')
    parser.add_argument("--name", action="store", help="name of sample. Used in naming files")
    parser.add_argument("-odir",  help="Output directory to save files to")
    args = parser.parse_args()
    return args

if __name__=="__main__":
    args = parse_args()
    annovar = pd.read_csv(args.annovar, sep = '\t')
    gnomad_cols = [c for c in annovar.columns if c.startswith('AF')]
    gnomad_map = {}
    for c in gnomad_cols:
        gnomad_map[c] = 'gnomAD_genomes_' + c
    annovar.rename(columns = gnomad_map, inplace=True)
    intervar = pd.read_csv(args.intervar, sep = '\t')
    snpeff = pd.read_csv(args.snpeff, dtype=object, sep = '\t')

    snpeff['ID'] = snpeff['ID'].apply(lambda x: x if x.startswith('rs') else '.')
    
    merged = merge_annotations(annovar, intervar, snpeff)
    outfile = os.path.join(args.odir, f"{args.name}.merged_annotations.tsv")
    merged.to_csv(f"{args.odir}/{args.name}_merged_annotations.tsv",  sep='\t', index=False)
    save_by_feature(merged, args.name, args.odir)