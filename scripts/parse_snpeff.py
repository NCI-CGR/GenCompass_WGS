#!/usr/bin/env python3

import argparse
from contextlib import contextmanager
from collections import defaultdict
import sys


class SnpEffANNSubField:
    """Parsed data of a single SnpEff ANN entry
    Details for each sub-field can be found in SnpEff documentation
    https://pcingola.github.io/SnpEff/se_inputoutput/ 
    """
    def __init__(self, annotation_subfield):
        self.raw_annotation = annotation_subfield.strip('ANN=')


        split_annotaion = self.raw_annotation.split('|')
        self.allele = split_annotaion[0].strip('ANN=')
        self.annotation = split_annotaion[1]
        self.putative_impact = split_annotaion[2]
        self.gene_name = split_annotaion[3]
        self.gene_id = split_annotaion[4]
        self.feature_type = split_annotaion[5]
        self.feature_id = split_annotaion[6]
        self.transcript_biotype = split_annotaion[7]
        self.rank__total = split_annotaion[8]
        self.hgvs_c = split_annotaion[9]
        self.hgvs_p = split_annotaion[10]
        self.cdna_position__cdna_length = split_annotaion[11]
        self.cds_position__cds_length = split_annotaion[12]
        self.protein_position__protein_length = split_annotaion[13]
        self.distance_to_feature = split_annotaion[14]
        self.errors = split_annotaion[15]

    def __repr__(self) -> str:
        return self.raw_annotation
    def __str__(self) -> str:
        return self.raw_annotation
    
    def is_high_impact(self):
        return self.putative_impact in ('MODERATE', 'HIGH')


class SnpEffField:
    def __init__(self, annotation):
        self.annotation = annotation.strip('ANN=')
        split_annotation = self.annotation.split(',')
        self.subfield_dict = defaultdict(list)
        for annot in split_annotation:
            subfield = SnpEffANNSubField(annot)
            self.subfield_dict[subfield.putative_impact].append(subfield)
        # self.subfields = [SnpEffANNSubField(x) for x in split_annotation]
        self.highest_impact = self.find_highest_impact()
    
    def __repr__(self) -> str:
        return self.annotation

    @property
    def gene(self):
        gene_set = set([subfield.gene_name for subfield in self.subfield_dict[self.highest_impact] if subfield.gene_name != ''])
        return ';'.join(gene_set)
    
    def find_highest_impact(self):
        impact_levels = {'MODIFIER': 0, 'LOW': 1, 'MODERATE': 2, 'HIGH': 3}
        current_impact_level = 'MODIFIER'
        for impact in self.subfield_dict.keys():
            if impact_levels[impact] > impact_levels[current_impact_level]:
                current_impact_level = impact
        return current_impact_level
    
    @property
    def hgvs_c(self):
        
        hgvs_c = [subfield.hgvs_c for subfield in self.subfield_dict[self.highest_impact] if subfield.hgvs_c != '']
        
        if len(hgvs_c) > 0:
            return ';'.join(hgvs_c)
        else:
            return '.'
    
    @property
    def hgvs_p(self):
        hgvs_p = [subfield.hgvs_p for subfield in self.subfield_dict[self.highest_impact] if subfield.hgvs_p != '']
        if len(hgvs_p) > 0:
            return ';'.join(hgvs_p)
        else:
            return '.'
    
    @property
    def effect(self):
        effect = [subfield.annotation for subfield in self.subfield_dict[self.highest_impact]]
        effect = list(set(effect))
        effect.sort()
        if len(effect) > 0:
            return ';'.join(effect)
        else:
            return '.'

        
    @property
    def num_impact_transcripts(self):
        return len(self.subfield_dict[self.highest_impact])
        # total = 0
        # for sub in self.subfields:
        #     if sub.is_high_impact():
        #         total = total + 1
        # return total
    
    @property
    def is_classic_high_impact(self):      
        classic_high_impact_classifications = ['frameshift_variant', 'stop_gained','stop_lost','start_lost','splice_acceptor_variant','splice_donor_variant']
        if self.highest_impact == 'HIGH' and any(x in self.effect for x in classic_high_impact_classifications):
            return 'Y'
        return 'N'
        

    @property
    def total_transcripts(self):
        total = 0
        for impact, transcripts in self.subfield_dict.items():
            total = total + len(transcripts)
        return total


@contextmanager
def file_or_stdout(file_name):
    if file_name is None:
        yield sys.stdout
    else:
        with open(file_name, 'w') as out_file:
            yield out_file

def parse_args():
    parser = argparse.ArgumentParser("Parse snpEff tsv.  Breaks up annotation to determine highest impact and the number of transcripts that have a MEDIUM or HIGH impact. Annotation and set information should already be extracted using bcftools query")
    parser.add_argument('--include-header', dest='include_header',action='store_true', help='Include header column names')
    parser.add_argument('--write-unannotated', dest='write_unannotated',action='store_true', help='Include variants where no snpEff annotation exists')
    parser.add_argument('infile', help='Input snpEff VCF file')
    parser.add_argument("outfile", nargs='?', help="Output snpEff tsv filename. If not provided then output to stdout")
    args = parser.parse_args()
    return args.infile, args.outfile, args.include_header, args.write_unannotated


def parse_snpEff(infile, outfile, include_header, write_unannotated):
    header = ''
    with open(infile) as ifile, file_or_stdout(outfile) as ofile:
        if include_header:
            ofile.write('CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tset\tANN\thighest_impact\tnum_impact_transcripts\ttotal_transcripts\tHGVS.c\tHGVS.p\teffect\tSnpEffGene\n')
        for line in ifile.readlines():
            if line.startswith('#'):
                # header = header + line
                continue

            chrom, pos, id, ref, alt, qual, filt, set, ann = line.strip('\n').split('\t')
            if ann != '.':
                
                ann =  SnpEffField(ann)
                oline = f'{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filt}\t{set}\t{ann}\t{ann.highest_impact}\t{ann.num_impact_transcripts}\t{ann.total_transcripts}\t{ann.hgvs_c}\t{ann.hgvs_p}\t{ann.effect}\t{ann.gene}\n'
                ofile.write(oline)
            elif write_unannotated:
                oline = f'{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filt}\t{set}' + ('\t.' * 8) + '\n'
                ofile.write(oline)



if __name__=="__main__":
     infile, outfile, include_header, write_unannotated = parse_args()
     parse_snpEff(infile, outfile, include_header, write_unannotated)

    


    
