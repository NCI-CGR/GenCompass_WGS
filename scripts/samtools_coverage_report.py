
import pandas as pd
import os
import click

CHROMOSOMES=[f'chr{x}' for x in range(1,23) ] + ['chrX', 'chrY', 'chrM']
CHROMOSOMES
# %%
def read_samtools_coverage_table(path):
    sample = os.path.basename(path).split('.')[0]
    print(sample)
    coverage = pd.read_csv(path, sep='\t')
    cols = ['Sample ID']+list(coverage.columns)
    coverage = coverage[coverage['#rname'].isin(CHROMOSOMES)]
    coverage['Sample ID'] = sample
    coverage = coverage[cols]
    return coverage


@click.command()
@click.option('--samtools_coverage_fof')
@click.option('--project_name')
@click.option('-odir', '--output_directory')
def samtools_coverage_report(samtools_coverage_fof, project_name, output_directory):
    os.makedirs(output_directory, exist_ok=True)
    output_path=f'{output_directory}/{project_name}.samtools_coverage_report.csv'
    if os.path.exists(output_path):
        os.remove(output_path)
    is_first_sample=True
    for path in open(samtools_coverage_fof):
        path = path.strip()
        coverage = read_samtools_coverage_table(path)
        coverage.to_csv(output_path, index=False, mode='a', header=is_first_sample)
        is_first_sample=False

    

if __name__=="__main__":
    samtools_coverage_report()