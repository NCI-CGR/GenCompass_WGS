#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import json
import random
import string
import datetime    
import re
import logging
from helpers import Manifest
import helpers
import pathlib

GENCOMPASS_PATH = pathlib.Path(__file__).parent.parent.resolve()

DEFAULT_OPTIONS_JSON={
        "delete_intermediate_output_files": True,
        "final_workflow_outputs_dir":"results",
        "use_relative_output_paths": True,
        "final_workflow_log_dir": "logging",
        "final_call_logs_dir": "logging"
    }
DEFAULT_RUN_PARAMETERS=['#SBATCH --cpus-per-task=8','#SBATCH --mem 25g', '#SBATCH --time 8:00:00', 'module load cromwell', 'module load singularity']
DEFAULT_SWARM_PARAMETERS=['#SWARM --threads-per-process 8', '#SWARM --gb-per-process 25',  '#SWARM --time 4:00:00', '#SWARM --module cromwell,singularity', '#SWARM --sbatch "--export SINGULARITY_CACHEDIR=/data/COVID_WGS/wgs_pipeline_performance/singularity_cache"']
PREMAP_QC = 'premap_qc'
MAPPING = 'mapping'
VARIANT_CALLING='variant_calling'
JOINT_GENOTYPE = 'joint_genotype'
HARMONIZE = 'harmonize'


def get_random_string(length: int) -> str:
    """
    Return a string of <length> random uppercase letters
    """
    return "".join(random.choice(string.ascii_uppercase) for i in range(0, length))

def get_random_run_id(length: int) -> str:
    today = datetime.datetime.now()
    run_id = f"{today.year}{today.month}{today.day}_{get_random_string(length)}"
    return run_id


def create_fastq_files_table(fastq_files_path: str, manifest: Manifest, samples = []) :
    def get_sample_id(fastq_path: str) -> str:
        fastq_filename = os.path.basename(fastq_path)
        sample_id  = fastq_filename.split('_')[0]
        return sample_id
    fastq_files  = pd.read_csv(fastq_files_path, header=None, names=['Fastq Location'])
    fastq_files['Sample Run ID'] = fastq_files['Fastq Location'].apply(get_sample_id)

    fastq_files['Sample ID'] = fastq_files['Sample Run ID'].map(manifest.get_sample_run_id_map())
    #   = fastq_files.merge(manifest[['Sample ID', 'Sample Run ID']], how ='inner', on = 'Sample Run ID')
    fastq_files  = fastq_files [['Sample ID', 'Fastq Location']]
    if samples is not None and len(samples) > 0:
        fastq_files = fastq_files[fastq_files['Sample ID'].isin(samples)]
    return fastq_files    


def get_sample_ids(args_sample_ids, sample_list, fastq_files):
    sample_ids = []
    if args_sample_ids is not None:
        sample_ids = sample_ids + args_sample_ids
    if sample_list is not None:
        sample_list = pd.read_csv(sample_list, sep='\t', header=None, names=['Sample ID'])
        sample_list = list(sample_list['Sample ID'].unique())
        sample_ids = sample_ids + sample_list
    if len(sample_ids) == 0:
        sample_ids = list(fastq_files['Sample ID'].unique())
    
    sample_ids = set(sample_ids)
    return sample_ids


class WorkflowInstance:
    def __init__(
        self, wdl_path: str, cromwell_invocation: str,
        sample_id, output_dir: str, input_map: dict, 
        options_map: dict,
        runtime_paramers: list):
        self.wdl_path: str = wdl_path
        self.runtime_paramers: list = runtime_paramers
        self.options_path: str = None
        self.input_path: str = None
        self.workflow_name: str = None
        self.cromwell_invocation: str = cromwell_invocation
        self.inputs_map: dict = input_map.copy()


        if options_map == {}:
            self.options_map = DEFAULT_OPTIONS_JSON.copy()
        else:
            self.options_map = options_map.copy()

        self.output_dir: str = output_dir

        if sample_id is None:
            self.set_random_sample_id()
        else:
            self.sample_id = sample_id

        self.set_auto_workflow_name()
        self.set_auto_run_path()
        self.set_auto_input_path()
        self.set_auto_options_path()
        self.set_cromwell_invocation()
    
    
    
    def set_random_sample_id(self):
        self.sample_id = get_random_string(8)

    def set_auto_workflow_name(self):
        workflow_name = os.path.basename(self.wdl_path)
        self.workflow_name = os.path.splitext(workflow_name)[0]
    
    def set_auto_run_path(self):
        if self.sample_id is None:
            self.set_random_sample_id()
        self.run_path = os.path.join(self.output_dir, f'{self.sample_id}.{self.workflow_name}_run.sh')

    def set_auto_input_path(self):
        if self.sample_id is None:
            self.set_random_sample_id()
        self.input_path = os.path.join(self.output_dir,f'{self.workflow_name}_inputs', f'{self.sample_id}.{self.workflow_name}_inputs.json')

    def set_auto_options_path(self):
        if self.sample_id is None:
            self.set_random_sample_id()
        self.options_path = os.path.join(self.output_dir,f'{self.workflow_name}_options', f'{self.sample_id}.{self.workflow_name}_options.json')

    
    def update_input_map(self, key, value):
        if value is not None:
            self.inputs_map[key] = value
        else:
             self.inputs_map.pop(key, None)
    
    def update_options_map(self, key, value):
        self.options_map[key] = value
    
    def set_cromwell_invocation(self):
        if self.sample_id is None:
            self.set_random_sample_id()
        if self.options_path is None:
            self.set_auto_options_path()
        if self.input_path is None:
            self.set_auto_input_path()
        
        
        self.cromwell_invocation = self.cromwell_invocation.replace("<<WORKFLOW>>", self.wdl_path)
        self.cromwell_invocation = self.cromwell_invocation.replace("<<INPUT_JSON>>", self.input_path)
        self.cromwell_invocation = self.cromwell_invocation.replace("<<OPTIONS_JSON>>", self.options_path)

    def save_inputs(self):
        json.dump(self.inputs_map, open(self.input_path, 'w'), sort_keys=True, indent=4)

    def save_options(self):
        json.dump(self.options_map, open(self.options_path, 'w'), sort_keys=True, indent=4)
    
    def save_run(self):
        with open(self.run_path, 'w') as run_writer:
            run_writer.write('#!/bin/bash\n')
            for line in self.runtime_paramers:
                run_writer.write(line + '\n')
            run_writer.write('\n')
            run_writer.write(self.cromwell_invocation)


class BatchBuilder:
    def __init__(
        self, project: str, run_id:str, wdl_path: str, 
        cromwell_invocation: str,
        input_template: dict, 
        output_dir:str,
          samples: set, options_template: dict, 
        runtime_parameters:list,
        environment:str):

        self.project: str = project
        self.run_id: str = run_id
        self.input_template: dict = input_template

        self.wdl_path: str = wdl_path
        self.cromwell_invocation=cromwell_invocation
        self.output_dir: str = output_dir
        self.samples: set = samples
        self.options_template: json = options_template

        self.runtime_parameters=runtime_parameters 
        self.environment=environment

        self.workflow_name: str = os.path.basename(self.wdl_path).strip('.wdl')

        inputs_dir=os.path.join(self.output_dir, f"{self.workflow_name}_inputs")
        options_dir=os.path.join(self.output_dir, f"{self.workflow_name}_options")
        os.makedirs(inputs_dir, exist_ok=True)
        os.makedirs(options_dir, exist_ok=True)

        if self.environment in ("swarm", "gcp", "aws"):
            ext = "swarm" if self.environment == "swarm" else "sh"
            self.batch_ofile = os.path.join(self.output_dir, f"{self.project}_{self.run_id}_{self.workflow_name}.{ext}")
            with open(self.batch_ofile, "w") as ofile:
                for line in self.runtime_parameters:
                    ofile.write(line + "\n")

    def create_workflow_files(self):
        pass
    
    def save_workflow(self, workflow:WorkflowInstance):
        if self.environment in ("slurm", "local"):
            workflow.save_run()
        elif self.environment in ("swarm", "gcp", "aws"):
            with open(self.batch_ofile, "a") as ofile:
                ofile.write(workflow.cromwell_invocation + "\n")


class PremapQCBatchBuilder(BatchBuilder):
    def __init__(self, project: str, run_id: str, wdl_path: str, cromwell_invocation: str, input_template: dict, output_dir: str, samples: set, options_template: dict, runtime_parameters: list, environment: str):
        super().__init__(project, run_id, wdl_path, cromwell_invocation, input_template, output_dir, samples, options_template, runtime_parameters, environment)

    
    def create_workflow_files(self):
        # self.fastq_files.to_csv(os.path.join(self.output_dir, self.fastq_table_filename), sep='\t', index=False)
        for sample_id in self.samples:
            sample_workflow = WorkflowInstance(wdl_path=self.wdl_path, cromwell_invocation=self.cromwell_invocation, sample_id=sample_id,
                                               output_dir=self.output_dir, input_map=self.input_template, 
                                               options_map=self.options_template, runtime_paramers=self.runtime_parameters)
            if self.environment == "aws":
                sample_workflow.update_input_map('sampleID', sample_id)
            else:
                sample_workflow.update_input_map('PremapQC.sampleID', sample_id)
            
            if self.environment in ['slurm', 'swarm', 'local']:
                if os.path.exists(os.path.join(self.options_template['final_workflow_outputs_dir'], f'premap_qc/fastp/{sample_id}')):
                    sample_workflow.update_input_map('PremapQC.runFastp', False)
                if os.path.exists(os.path.join(self.options_template['final_workflow_outputs_dir'], f'premap_qc/fastqc/{sample_id}')):
                    sample_workflow.update_input_map('PremapQC.runFastQC', False)
                if os.path.exists(os.path.join(self.options_template['final_workflow_outputs_dir'], f'premap_qc/fastq_screen/{sample_id}')):
                    sample_workflow.update_input_map('PremapQC.runFastqScreen', False)
                    

            sample_workflow.update_options_map('final_workflow_log_dir', os.path.join(sample_workflow.options_map['final_workflow_log_dir'], f'premap-qc/{sample_id}') )
            sample_workflow.update_options_map('final_call_logs_dir', os.path.join(sample_workflow.options_map['final_call_logs_dir'], f'premap-qc/{sample_id}') )
            sample_workflow.save_inputs()
            sample_workflow.save_options()
            self.save_workflow(sample_workflow)
             
           
class MappingBatchBuilder(BatchBuilder):
    def __init__(self, project: str, run_id: str, wdl_path: str, cromwell_invocation: str, input_template: dict, output_dir: str, samples: set, options_template: dict, runtime_parameters: list, environment: str, fastq_files: pd.DataFrame):
        super().__init__(project, run_id, wdl_path, cromwell_invocation, input_template, output_dir, samples, options_template, runtime_parameters, environment)
        self.fastq_files = fastq_files
        
    
    def create_fastp_fastq_results_files(self):
        def get_fastp_result_loc(sample: str, original_fastq_filename: str):
            lane_id = helpers.extract_sample_lane_id(original_fastq_filename)
            paired_end = helpers.extract_sample_paired_end(original_fastq_filename)
            fastp_filename = os.path.join(self.options_template['final_workflow_outputs_dir'], f'premap_qc/fastp/{sample}/{lane_id}_{paired_end}_fastp.fastq.gz')
            return fastp_filename
            
        for sample_id in self.samples:
            sample_fastq = self.fastq_files[self.fastq_files['Sample ID'] == sample_id].copy()
            sample_fastq['fastp_filename'] = sample_fastq.apply(lambda row: get_fastp_result_loc(row['Sample ID'], row['Fastq Location']), axis = 1)
            sample_fastq[['fastp_filename']].to_csv(os.path.join(self.output_dir, 'mapping_inputs',f'{sample_id}.fastp_fastq_files.txt'), header=False, index=False, sep='\t')

    def create_workflow_files(self):
        
        def get_fastp_result_loc(sample: str, original_fastq_filename: str):
            lane_id = helpers.extract_sample_lane_id(original_fastq_filename)
            paired_end = helpers.extract_sample_paired_end(original_fastq_filename)
            fastp_filename = os.path.join(self.options_template['final_workflow_outputs_dir'], f'premap_qc/fastp/{sample}/{lane_id}_{paired_end}_fastp.fastq.gz')
            return fastp_filename
        for sample_id in self.samples:

            sample_workflow = WorkflowInstance(
                wdl_path=self.wdl_path,
                cromwell_invocation=self.cromwell_invocation,
                sample_id=sample_id,
                output_dir=self.output_dir,
                input_map=self.input_template,
                options_map=self.options_template,
                runtime_paramers=self.runtime_parameters
            )
            # sample_workflow.update_input_map('Mapping.sampleFastqLocations', os.path.join(self.output_dir, 'mapping_inputs',f'{sample_id}.fastp_fastq_files.txt'))
            sample_fastq = self.fastq_files[self.fastq_files['Sample ID'] == sample_id]
            sample_file_locations = [get_fastp_result_loc(row['Sample ID'], row['Fastq Location']) for _, row in sample_fastq.iterrows()]
            if self.environment == "aws":
                sample_workflow.update_input_map('sampleFastqFiles', sample_file_locations)
                sample_workflow.update_input_map('sampleID', sample_id)
            else:
                sample_workflow.update_input_map('Mapping.sampleFastqFiles', sample_file_locations)
                sample_workflow.update_input_map('Mapping.mappedBAM', None)

                if self.environment in ['slurm', 'swarm', 'local']:
                    if os.path.exists(os.path.join(self.options_template['final_workflow_outputs_dir'],f'fq2bam/{sample_id}/{sample_id}.bam')):
                        sample_workflow.update_input_map('Mapping.mappedBAM', os.path.join(self.options_template['final_workflow_outputs_dir'],f'fq2bam/{sample_id}/{sample_id}.bam'))
                        sample_workflow.update_input_map('Mapping.sampleFastqFiles', None)
                    
                    # bammetrics check
                    if os.path.exists(os.path.join(self.options_template['final_workflow_outputs_dir'],f'mapping_qc/bammetrics/{sample_id}.bammetrics.txt')):
                        sample_workflow.update_input_map('Mapping.runBamMetrics', False)

                    # collectmultiplemetrics check
                    if os.path.exists(os.path.join(self.options_template['final_workflow_outputs_dir'],f'mapping_qc/collectmultiplemetrics/{sample_id}')):
                        sample_workflow.update_input_map('Mapping.runCollectMultipleMetrics', False)

                    # kraken check
                    if os.path.exists(os.path.join(self.options_template['final_workflow_outputs_dir'],f'mapping_qc/kraken2/{sample_id}')):
                        sample_workflow.update_input_map('Mapping.runKraken', False)

                    # samtools coverage check
                    if os.path.exists(os.path.join(self.options_template['final_workflow_outputs_dir'],f'mapping_qc/samtools_coverage/{sample_id}.samtools_coverage.txt')):
                        sample_workflow.update_input_map('Mapping.runSamtoolsCoverage', False)
                
                    # somalier check
                    if os.path.exists(os.path.join(self.options_template['final_workflow_outputs_dir'],f'mapping_qc/somalier/extract/{sample_id}.somalier')):
                        sample_workflow.update_input_map('Mapping.runSomalier', False)

                    # verifybamid check
                    if os.path.exists(os.path.join(self.options_template['final_workflow_outputs_dir'],f'mapping_qc/verifybamid/{sample_id}')):
                        sample_workflow.update_input_map('Mapping.runVerifyBamID', False)
                
                sample_workflow.update_input_map('Mapping.sampleID', sample_id)

            sample_workflow.update_options_map('final_workflow_log_dir', os.path.join(sample_workflow.options_map['final_workflow_log_dir'],f'mapping/{sample_id}') )
            sample_workflow.update_options_map('final_call_logs_dir', os.path.join(sample_workflow.options_map['final_call_logs_dir'], f'mapping/{sample_id}') )
            sample_workflow.save_inputs()
            sample_workflow.save_options()
            self.save_workflow(sample_workflow)

class VariantCallingBatchBuilder(BatchBuilder):
    def __init__(self, project: str, run_id: str, wdl_path: str, cromwell_invocation: str, input_template: dict, output_dir: str, samples: set, options_template: dict, runtime_parameters: list, environment: str):
        super().__init__(project, run_id, wdl_path, cromwell_invocation, input_template, output_dir, samples, options_template, runtime_parameters, environment)
    
    def create_workflow_files(self):
        for sample_id in self.samples:
            sample_workflow = WorkflowInstance(
                wdl_path=self.wdl_path, 
                cromwell_invocation=self.cromwell_invocation,
                sample_id=sample_id, 
                output_dir=self.output_dir, 
                input_map=self.input_template, 
                options_map=self.options_template,
                runtime_paramers=self.runtime_parameters)
            if self.environment == "aws":
                sample_workflow.update_input_map('sampleID', sample_id)
                sample_workflow.update_input_map('sampleBAM', os.path.join(self.options_template['final_workflow_outputs_dir'], f'fq2bam/{sample_id}/{sample_id}.bam'))
                sample_workflow.update_input_map('sampleBAI', os.path.join(self.options_template['final_workflow_outputs_dir'], f'fq2bam/{sample_id}/{sample_id}.bam.bai'))
                sample_workflow.update_input_map('sampleBQSR', os.path.join(self.options_template['final_workflow_outputs_dir'], f'fq2bam/{sample_id}/{sample_id}.BQSR-REPORT.txt'))

            else:
                sample_workflow.update_input_map('VariantCalling.sampleID', sample_id)
                sample_workflow.update_input_map('VariantCalling.sampleBAM', os.path.join(self.options_template['final_workflow_outputs_dir'], f'fq2bam/{sample_id}/{sample_id}.bam'))
                sample_workflow.update_input_map('VariantCalling.sampleBAI', os.path.join(self.options_template['final_workflow_outputs_dir'], f'fq2bam/{sample_id}/{sample_id}.bam.bai'))
                sample_workflow.update_input_map('VariantCalling.sampleBQSR', os.path.join(self.options_template['final_workflow_outputs_dir'], f'fq2bam/{sample_id}/{sample_id}.BQSR-REPORT.txt'))

            sample_workflow.update_options_map('final_workflow_log_dir', os.path.join(sample_workflow.options_map['final_workflow_log_dir'],f'variant-calling/{sample_id}') )
            sample_workflow.update_options_map('final_call_logs_dir', os.path.join(sample_workflow.options_map['final_call_logs_dir'], f'variant-calling/{sample_id}') )
            sample_workflow.save_inputs()
            sample_workflow.save_options()
            self.save_workflow(sample_workflow)
        

class JointGenotypeBatchBuilder(BatchBuilder):
    def __init__(self, project: str, run_id: str, wdl_path: str, cromwell_invocation: str, input_template: dict, output_dir: str, samples: set, options_template: dict, runtime_parameters: list, environment: str, strelka_glnexus_config: str):
        super().__init__(project, run_id, wdl_path, cromwell_invocation, input_template, output_dir, samples, options_template, runtime_parameters, environment)
        self.callers = ['deepvariant', 'haplotypecaller', 'strelka2']
        self.strelka_glnexus_config = strelka_glnexus_config

    def create_gvcf_fof(self):
        for caller in self.callers:
            filename = os.path.join(self.output_dir,f"{self.workflow_name}_inputs", f'{self.project}.{self.run_id}.{caller}_vcf_files.txt')
            with open(filename, 'w') as ofile:
                for sample_id in self.samples:
                    if caller == 'strelka2':
                        path = os.path.join(self.options_template['final_workflow_outputs_dir'], f'{caller}/{sample_id}/{sample_id}.strelka.genome.vcf.gz')
                    else:
                        path = os.path.join(self.options_template['final_workflow_outputs_dir'], f'{caller}/{sample_id}/{sample_id}.{caller}.g.vcf.gz')
                    ofile.write(path + '\n')
    
    def create_workflow_files(self, mode="run"):
        for caller in self.callers:
            sample_workflow = WorkflowInstance(
                wdl_path=self.wdl_path, 
                cromwell_invocation=self.cromwell_invocation,
                sample_id=caller, 
                output_dir=self.output_dir, 
                input_map=self.input_template, 
                options_map=self.options_template,
                runtime_paramers=self.runtime_parameters)
            

            sample_workflow.update_input_map('JointGenotype.caller', caller)
            input_files_filename = os.path.join(self.output_dir,f"{self.workflow_name}_inputs", f'{self.project}.{self.run_id}.{caller}_vcf_files.txt')

            glnexus_config = {'haplotypecaller': 'gatk',
                              'deepvariant': 'DeepVariant',
                              'strelka2': self.strelka_glnexus_config }
            
            sample_workflow.update_input_map('JointGenotype.glnexusConfig', glnexus_config[caller])
            
            if mode =="run":
                sample_workflow.update_input_map('JointGenotype.variantFOF', input_files_filename)
                sample_workflow.inputs_map.pop('JointGenotype.calledVariants', None)
            
            elif mode =="server":
                with open(input_files_filename) as ifile:
                    gvcf_files = ifile.read().splitlines()
                    sample_workflow.update_input_map('JointGenotype.calledVariants', gvcf_files)
                    sample_workflow.inputs_map.pop('JointGenotype.variantFOF', None)
                os.remove(input_files_filename)
            elif mode =="aws-omics":
                with open(input_files_filename) as ifile:
                    gvcf_files = ifile.read().splitlines()
                    sample_workflow.update_input_map('calledVariants', gvcf_files)
                    sample_workflow.inputs_map.pop('variantFOF', None)
                os.remove(input_files_filename)




            sample_workflow.update_options_map('final_workflow_log_dir', os.path.join(sample_workflow.options_map['final_workflow_log_dir'],f'joint-genotype/{caller}') )
            sample_workflow.update_options_map('final_call_logs_dir', os.path.join(sample_workflow.options_map['final_call_logs_dir'], f'joint-genotype/{caller}') )
            sample_workflow.save_inputs()
            sample_workflow.save_options()
            self.save_workflow(sample_workflow)


class HarmonizeBatchBuilder(BatchBuilder):
    def __init__(self, project: str, run_id: str, wdl_path: str, cromwell_invocation: str, input_template: dict, output_dir: str, samples: set, options_template: dict, runtime_parameters: list, environment: str, bin_interval_names_file: str= None):
        super().__init__(project, run_id, wdl_path, cromwell_invocation, input_template, output_dir, samples, options_template, runtime_parameters, environment)
        self.bin_interval_names_file = bin_interval_names_file
    def create_workflow_files(self):
        if self.bin_interval_names_file is None:
            self.create_workflow_files_concat()
        else:
            self.create_workflow_files_intervals()

    def create_workflow_files_concat(self):
        sample_workflow = WorkflowInstance(
                wdl_path=self.wdl_path, 
                cromwell_invocation=self.cromwell_invocation,
                sample_id=f'{self.project}_{self.run_id}', 
                output_dir=self.output_dir, 
                input_map=self.input_template, 
                options_map=self.options_template,
                runtime_paramers=self.runtime_parameters)
        if self.environment == "aws":
            sample_workflow.update_input_map('haplotypecallerVCF', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/haplotypecaller/haplotypecaller.concat.bcf.gz'))
            sample_workflow.update_input_map('haplotypecallerVCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/haplotypecaller/haplotypecaller.concat.bcf.gz.tbi'))

            sample_workflow.update_input_map('deepvariantVCF', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/deepvariant/deepvariant.concat.bcf.gz'))
            sample_workflow.update_input_map('deepvariantVCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/deepvariant/deepvariant.concat.bcf.gz.tbi'))

            sample_workflow.update_input_map('strelka2VCF', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/strelka2/strelka2.concat.bcf.gz'))
            sample_workflow.update_input_map('strelka2VCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/strelka2/strelka2.concat.bcf.gz.tbi'))
        else:
            sample_workflow.update_input_map('haplotypecallerVCF', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/haplotypecaller/haplotypecaller.concat.bcf.gz'))
            sample_workflow.update_input_map('haplotypecallerVCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/haplotypecaller/haplotypecaller.concat.bcf.gz.tbi'))

            sample_workflow.update_input_map('deepvariantVCF', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/deepvariant/deepvariant.concat.bcf.gz'))
            sample_workflow.update_input_map('deepvariantVCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/deepvariant/deepvariant.concat.bcf.gz.tbi'))

            sample_workflow.update_input_map('strelka2VCF', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/strelka2/strelka2.concat.bcf.gz'))
            sample_workflow.update_input_map('strelka2VCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], 'joint_genotype/strelka2/strelka2.concat.bcf.gz.tbi'))
            

        sample_workflow.update_options_map('final_workflow_log_dir', os.path.join(sample_workflow.options_map['final_workflow_log_dir'],f'harmonize') )
        sample_workflow.update_options_map('final_call_logs_dir', os.path.join(sample_workflow.options_map['final_call_logs_dir'], f'harmonize') )
        sample_workflow.save_inputs()
        sample_workflow.save_options()
        self.save_workflow(sample_workflow)
   
    def create_workflow_files_intervals(self):
        for interval in open(self.bin_interval_names_file):
            interval = interval.strip()
            sample_workflow = WorkflowInstance(
                    wdl_path=self.wdl_path, 
                    cromwell_invocation=self.cromwell_invocation,
                    sample_id=f'{self.project}_{self.run_id}_{interval}', 
                    output_dir=self.output_dir, 
                    input_map=self.input_template, 
                    options_map=self.options_template,
                    runtime_paramers=self.runtime_parameters)
            
            if self.environment == "aws":
                sample_workflow.update_input_map('project', f'{self.input_template["project"]}_{interval}')
                sample_workflow.update_input_map('haplotypecallerVCF', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/haplotypecaller/haplotypecaller.{interval}.bcf.gz'))
                sample_workflow.update_input_map('haplotypecallerVCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/haplotypecaller/haplotypecaller.{interval}.bcf.gz.tbi'))

                sample_workflow.update_input_map('deepvariantVCF', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/deepvariant/deepvariant.{interval}.bcf.gz'))
                sample_workflow.update_input_map('deepvariantVCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/deepvariant/deepvariant.{interval}.bcf.gz.tbi'))

                sample_workflow.update_input_map('strelka2VCF', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/strelka2/strelka2.{interval}.bcf.gz'))
                sample_workflow.update_input_map('strelka2VCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/strelka2/strelka2.{interval}.bcf.gz.tbi'))
            
            else:
                sample_workflow.update_input_map('Harmonize.project', f'{self.input_template["Harmonize.project"]}_{interval}')
                sample_workflow.update_input_map('Harmonize.haplotypecallerVCF', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/haplotypecaller/haplotypecaller.{interval}.bcf.gz'))
                sample_workflow.update_input_map('Harmonize.haplotypecallerVCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/haplotypecaller/haplotypecaller.{interval}.bcf.gz.tbi'))

                sample_workflow.update_input_map('Harmonize.deepvariantVCF', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/deepvariant/deepvariant.{interval}.bcf.gz'))
                sample_workflow.update_input_map('Harmonize.deepvariantVCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/deepvariant/deepvariant.{interval}.bcf.gz.tbi'))

                sample_workflow.update_input_map('Harmonize.strelka2VCF', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/strelka2/strelka2.{interval}.bcf.gz'))
                sample_workflow.update_input_map('Harmonize.strelka2VCFIndex', os.path.join(self.options_template['final_workflow_outputs_dir'], f'joint_genotype/strelka2/strelka2.{interval}.bcf.gz.tbi'))
                
            sample_workflow.update_options_map('final_workflow_log_dir', os.path.join(sample_workflow.options_map['final_workflow_log_dir'],f'harmonize/{interval}') )
            sample_workflow.update_options_map('final_call_logs_dir', os.path.join(sample_workflow.options_map['final_call_logs_dir'], f'harmonize/{interval}') )
            sample_workflow.save_inputs()
            sample_workflow.save_options()
            self.save_workflow(sample_workflow)
   

def set_options(options_file):
    if options_file:
        options = json.load( open(options_file))
    else:
        options = DEFAULT_OPTIONS_JSON
    
    if 'final_workflow_outputs_dir' not in options.keys():
        options['final_workflow_outputs_dir'] = DEFAULT_OPTIONS_JSON['final_workflow_outputs_dir']
    if 'final_workflow_log_dir' not in options.keys():
        options['final_workflow_log_dir'] = DEFAULT_OPTIONS_JSON['final_workflow_log_dir']
    if 'final_call_logs_dir' not in options.keys():
        options['final_call_logs_dir'] = DEFAULT_OPTIONS_JSON['final_call_logs_dir']
    return options


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--project', action='store', required=True, dest='project', help='Name of project. Used in file naming')
    parser.add_argument('--run_id', action='store', required=False, dest='run_id', help='Run ID of subsample. Used in file naming. If not provided then random run_id assigned')
    parser.add_argument("--cromwell_invocation", action="store")
    parser.add_argument('--workflows', 
        choices=[PREMAP_QC, MAPPING, VARIANT_CALLING, JOINT_GENOTYPE, HARMONIZE], 
        default=[PREMAP_QC, MAPPING, VARIANT_CALLING, JOINT_GENOTYPE, HARMONIZE], 
        nargs='+', help='Workflows to run')
    parser.add_argument('--fastq_list', action='store', dest='fastq_list', required=True, help='Fastq file of files. List of files separated by lines')
    parser.add_argument('--manifest', required=True, help='Manifest file. Contains Sample ID and Flowcell. Optionally includes Sample Run ID')

    parser.add_argument('--sample_ids', required=False, nargs='*', help='Sample IDs. If sample IDs are provided only inputs and run scripts for those samples will be generated')
    parser.add_argument('--sample_list', required=False, action='store', help='File of sample ids to process. If sample_ids and sample_list blank then all samples in fastq file will be generated.')


    parser.add_argument('--premap_qc_inputs_template',  action = "store", required=False, help='PremapQC inputs template file. JSON file')
    parser.add_argument('--mapping_inputs_template', action = "store", required=False, help='Mapping inputs template file. JSON file')
    parser.add_argument('--variant_calling_inputs_template', action='store', required=False, help='Variant calling inputs template file. JSON file')
    parser.add_argument('--joint_genotype_inputs_template', action='store', required=False, help='Joint genotype inputs template file. JSON file')
    parser.add_argument('--harmonize_inputs_template', action='store', required=False, help='Harmonize inputs template file. JSON file')
    parser.add_argument('--bin_interval_names', action='store', required=False, help='Bin interval names - used in naming joint genotype and harmonize outputs')
    parser.add_argument('--joint_genotype_strelka_glnexus_config', action='store', required=False, help='Strelka glnexus config file')

    parser.add_argument('--options_json', dest='options_json', required=False, help='Options file to use. If not provided a default options JSON is used.')

    parser.add_argument('--runtime_parameters', nargs='*', default=[], help='Runtime parameters for slurm or swarm file. If in the biowulf ecosystem it is necessary to load modules and set the singularity cache directory')
    parser.add_argument('--environment', choices=['slurm', 'swarm', 'local', 'aws', 'gcp'], default='local', help='System and environment to use. Slurm generates individual run.sh files while swarm generates one swarm file per workflow')
    parser.add_argument("--mode", choices=['run', 'server'], default='run')

    
    parser.add_argument('-odir', '--output_directory', dest='odir', default ='./', help='Output directory')

    args = parser.parse_args()
    return args


def get_inputs(workflow, input_template_path, environment):
    if input_template_path is None:
            inputs_template=f"input_templates/{workflow}.inputs_template.json"
            logging.info(f'{workflow}input not provided. Using {inputs_template}')
    else:
        inputs_template = input_template_path
    inputs_template = json.load(open(inputs_template))
    runtime_folder = 'hpc' if environment in ['slurm', 'swarm', 'local'] else environment
    runtime_inputs = json.load(open(os.path.join(GENCOMPASS_PATH, f'templates/{runtime_folder}/{workflow}.runtime_inputs.json')))
    inputs_template = inputs_template | runtime_inputs
    return inputs_template

def main():
    args = parse_args()
    
    

    if args.run_id is None:
        run_id = get_random_run_id(8)
    else:
        run_id = args.run_id
    
    options_template = set_options(args.options_json)
    
    

    # manifest = read_manifest(args.manifest)
    manifest = Manifest(args.manifest)
    fastq_files = create_fastq_files_table(args.fastq_list, manifest)

    sample_ids = get_sample_ids(args.sample_ids, args.sample_list, fastq_files)


    if PREMAP_QC in args.workflows:
        logging.info('Creating premap-qc files')
        premap_odir = os.path.join(args.odir, 'premap-qc')
        premap_qc_wdl = os.path.join(GENCOMPASS_PATH, 'workflows/premap_qc.wdl')

        premap_qc_inputs_template = get_inputs('premap_qc', args.premap_qc_inputs_template, args.environment)
        
        premap_builder = PremapQCBatchBuilder(
            project=args.project, run_id=run_id, samples=sample_ids,
            wdl_path=premap_qc_wdl,
            cromwell_invocation=args.cromwell_invocation,
            input_template=premap_qc_inputs_template, options_template=options_template,
            runtime_parameters=args.runtime_parameters, environment=args.environment,
            output_dir=args.odir
        )
        premap_builder.create_workflow_files()
        
    if MAPPING in args.workflows:
        logging.info('Creating mapping files')
        mapping_odir = os.path.join(args.odir, 'mapping')
        mapping_wdl = os.path.join(GENCOMPASS_PATH, 'workflows/mapping.wdl')

        mapping_inputs_template = get_inputs('mapping', args.mapping_inputs_template, args.environment)

        mapping_builder = MappingBatchBuilder(
            project=args.project, run_id=run_id, wdl_path=mapping_wdl,
            cromwell_invocation=args.cromwell_invocation,
            samples=sample_ids,
            input_template=mapping_inputs_template, 
            options_template = options_template,
            runtime_parameters=args.runtime_parameters,
            environment=args.environment, 
            fastq_files=fastq_files,
            output_dir=args.odir
        )
        # mapping_builder.create_fastp_fastq_results_files()
        mapping_builder.create_workflow_files()
    
    if VARIANT_CALLING in args.workflows:
        logging.info('Creating variant-calling run files')
        variant_calling_odir = os.path.join(args.odir, 'variant-calling')
        variant_calling_wdl = os.path.join(GENCOMPASS_PATH, 'workflows/variant_calling.wdl')
 
        variant_calling_inputs_template = get_inputs('variant_calling', args.variant_calling_inputs_template, args.environment)
        variant_calling_builder = VariantCallingBatchBuilder(
            project=args.project, run_id=run_id, wdl_path=variant_calling_wdl,
            cromwell_invocation=args.cromwell_invocation,
            samples=sample_ids,
            input_template=variant_calling_inputs_template, 
            options_template = options_template, 
            output_dir=args.odir,
            runtime_parameters=args.runtime_parameters,
            environment=args.environment
            )
        variant_calling_builder.create_workflow_files()
    
    if JOINT_GENOTYPE in args.workflows:
        logging.info('Creating joint genotype run files')
        

        joint_genotype_odir= os.path.join(args.odir, 'joint-genotype')
        joint_genotype_wdl = os.path.join(GENCOMPASS_PATH, 'workflows/joint_genotype.wdl')

        joint_genotype_inputs_template = get_inputs('joint_genotype', args.joint_genotype_inputs_template, args.environment)

        
        if args.joint_genotype_strelka_glnexus_config is None:
            strelka_glnexus_config = os.path.join(GENCOMPASS_PATH, 'config/strelka2_glnexus.yml')
        else:
            strelka_glnexus_config = args.joint_genotype_strelka_glnexus_config

        joint_genotype_builder = JointGenotypeBatchBuilder(
            project=args.project, run_id=run_id, wdl_path=joint_genotype_wdl,
            cromwell_invocation=args.cromwell_invocation,
            samples=sample_ids,
            input_template=joint_genotype_inputs_template, options_template = options_template, 
            output_dir=args.odir,
            runtime_parameters=args.runtime_parameters,
            strelka_glnexus_config=strelka_glnexus_config,
            environment=args.environment
            )
        
        joint_genotype_builder.create_gvcf_fof()
        joint_genotype_builder.create_workflow_files(mode=args.mode)


    if HARMONIZE in args.workflows:
        logging.info('Creating harmonize run files')
        
        harmonize_odir= os.path.join(args.odir, 'harmonize')
        harmonize_wdl = os.path.join(GENCOMPASS_PATH, 'workflows/harmonize.wdl')

        harmonize_inputs_template = get_inputs('harmonize', args.harmonize_inputs_template, args.environment)

        harmonize_builder = HarmonizeBatchBuilder(
            project=args.project, run_id=run_id, wdl_path=harmonize_wdl,
            cromwell_invocation=args.cromwell_invocation,
            samples=sample_ids,
            input_template=harmonize_inputs_template, options_template = options_template, 
            output_dir=args.odir,
            runtime_parameters=args.runtime_parameters,
            environment=args.environment,
            bin_interval_names_file=args.bin_interval_names
            )
        harmonize_builder.create_workflow_files()

if __name__ =="__main__":
    main()