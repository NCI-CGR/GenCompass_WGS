#!/bin/bash
ANNOVAR_HOME=Software/annovar



# Convert VCF to AVINPUT
# Converting allows for output to be csv format, which allows for easier merging of samples 
BUILDVER=hg38
VCF=covnet_500sample.left.chr22.vcf
AVINPUT=covnet_500sample.left.chr22.avinput


perl $ANNOVAR_HOME/convert2annovar.pl -format vcf4 --keepindelref $VCF > $AVINPUT


# Run ANNOVAR
# Output is file annovar.chr22.hg19_multianno.txt
perl $ANNOVAR_HOME/table_annovar.pl $AVINPUT $ANNOVAR_DATA/$BUILDVER --buildver $BUILDVER --out annovar.chr22 --remove --protocol refGene,gnomad30_genome,clinvar_20220320,cytoBand --operation g,f,f,r --nastring .


# Run INTERVAR
InterVar \
-i $AVINPUT \
-d $ANNOVAR_DATA/$BUILDVER \
-b $BUILDVER \
-o intervar_chr22_avinput

# Run SnpSift on gnomad - not implemented yet
GNOMAD_DB="path/to/gnomad"
java -jar $SNPSIFT_JAR annotate -info AF,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth -name gnomAD_genomes_  $GNOMAD_DB $VCF > covnet_chr22.gnomad.vcf

# Run snpEff
java -jar $SNPEFF_JAR \
-c $SNPEFF_HOME/snpEff.config \
-v -no-intergenic \
-no INTRAGENIC \
-no-downstream \
-no-upstream \
-no-utr \
-no PROTEIN_STRUCTURAL_INTERACTION_LOCUS \
-no NEXT_PROT \
-lof -noStats \
GRCh38.86 \
$VCF > covnet_chr22.snpEff.vcf

# Add dbSNP id field
DBSNP=/fdb/GATK_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz
# DBSNP=/fdb/genomebrowser/gbdb/hg38/snp/dbSnp153.bb
java -jar $SNPSIFT_JAR annotate -id $DBSNP covnet_chr22.snpEff.vcf > covnet_chr22.snpEff.dbSNPID.vcf

# Extract ANN from snpEff results
bcftools query  --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/set\t%INFO/ANN\n' \
covnet_chr22.snpEff.dbSNPID.vcf > covnet_chr22.snpEff.tsv

#-------------------------------------------
# Custom scripts for extracting relevant data 
# and merging data into single file
#-------------------------------------------
GENCOMPASS=/home/jordanbt/GenCompass

# Parse InterVar
python $GENCOMPASS/scripts/parse_intervar.py \
--intervar intervar_chr22_avinput.hg38_multianno.txt.intervar \
--output covnet_chr22.intervar_parsed.tsv

# Parse ANN field from snpEff and extract information
python $GENCOMPASS/scripts/parse_snpeff.py \
--include-header \
--write-unannotated \
covnet_chr22.snpEff.tsv  \
covnet_chr22.snpEff_parsed.tsv

# Merge samples
mkdir merged
python $GENCOMPASS/scripts/merge_annotations.py \
--annovar annovar.chr22.hg38_multianno.txt \
--intervar covnet_chr22.intervar_parsed.tsv \
--snpeff covnet_chr22.snpEff_parsed.tsv \
--name covnet_chr22 \
-odir merged

