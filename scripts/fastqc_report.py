#!/usr/bin/env python3

import pandas as pd
import click
import json
from helpers import Manifest
import statistics
import os



@click.command()
@click.option('--premap_multiqc_json')
@click.option('--manifest')
@click.option('--project_name')
@click.option('-odir', '--output_directory')
def fastqc_report(premap_multiqc_json, manifest, project_name, output_directory):
    with open(premap_multiqc_json, "r") as multiqc_json:
        json_data = json.load(multiqc_json)
    fastqc_data = json_data['report_saved_raw_data']['multiqc_fastqc']
    fastqc_lane = pd.DataFrame(fastqc_data).T.reset_index().rename(columns={'index': 'Lane ID'})
    columns_ordered = ['Sample ID','Sample Run ID'] + list(fastqc_lane.columns)
    fastqc_lane['Sample Run ID'] = fastqc_lane['Lane ID'].apply(lambda x: x.split('_')[0])
    

    sample_run_id_map = Manifest(manifest).get_sample_run_id_map()
    
    fastqc_lane['Sample ID'] = fastqc_lane['Sample Run ID'].map(sample_run_id_map)

    fastqc_lane = fastqc_lane[columns_ordered]

    def worst_result(results):
        if 'fail' in results.values:
            return 'fail'
        elif 'warn' in results.values:
            return 'warn'
        else:
            return 'pass'
    
    def mode_result(results):
        preference_order = ['pass', 'warn', 'fail']
        mode = statistics.multimode(results)
        sorted_mode = sorted(mode, key=lambda x: preference_order.index(x))
        return sorted_mode[0]

        


    fastqc_sample = fastqc_lane.groupby('Sample ID').agg(
        total_sequences=('Total Sequences', 'sum'),
        sequences_flagged_as_poor_quality=('Sequences flagged as poor quality', 'sum'),
        average_percent_gc=('%GC', 'mean'),
        average_total_deduplicated_percentage=('total_deduplicated_percentage', 'mean'),

        mode_base_statistics=('basic_statistics', mode_result),
        mode_per_base_sequence_quality=('per_base_sequence_quality', mode_result),
        mode_per_tile_sequence_quality=('per_tile_sequence_quality', mode_result),
        mode_per_sequence_quality_scores=('per_sequence_quality_scores', mode_result),
        mode_per_base_sequence_content=('per_base_sequence_content', mode_result),
        mode_per_sequence_gc_content=('per_sequence_gc_content', mode_result),
        mode_per_base_n_content=('per_base_n_content', mode_result),
        mode_sequence_length_distribution=('sequence_length_distribution', mode_result),
        mode_sequence_duplication_levels=('sequence_duplication_levels', mode_result),
        mode_overrepresented_sequences=('overrepresented_sequences', mode_result),
        mode_adapter_content=('adapter_content', mode_result),

        worst_base_statistics=('basic_statistics', worst_result),
        worst_per_base_sequence_quality=('per_base_sequence_quality', worst_result),
        worst_per_tile_sequence_quality=('per_tile_sequence_quality', worst_result),
        worst_per_sequence_quality_scores=('per_sequence_quality_scores', worst_result),
        worst_per_base_sequence_content=('per_base_sequence_content', worst_result),
        worst_per_sequence_gc_content=('per_sequence_gc_content', worst_result),
        worst_per_base_n_content=('per_base_n_content', worst_result),
        worst_sequence_length_distribution=('sequence_length_distribution', worst_result),
        worst_sequence_duplication_levels=('sequence_duplication_levels', worst_result),
        worst_overrepresented_sequences=('overrepresented_sequences', worst_result),
        worst_adapter_content=('adapter_content', worst_result),
        
    )
    fastqc_sample.reset_index(inplace=True)
    os.makedirs(output_directory, exist_ok=True)
    with pd.ExcelWriter(f"{output_directory}/{project_name}.pre-mapping-QC-C-fastqc_report.xlsx", mode="w", engine="openpyxl") as writer:
     
        fastqc_sample.to_excel(writer, sheet_name="FastQC Sample Report", index=False)
        fastqc_lane.to_excel(writer, sheet_name="FastQC Lane Report", index=False)
    return fastqc_sample



if __name__=="__main__":
    fastqc_report()