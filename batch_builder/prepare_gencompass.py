import click
import os
import shutil
import json
import pandas as pd
from helpers import Manifest
from bin_intervals import get_interval_names
import helpers

# Usage: prepare_gencompass.py [OPTIONS] COMMAND1 [ARGS]... [COMMAND2
#                              [ARGS]...]...

# Options:
#   --help  Show this message and exit.

# Commands:
#   all              Premap, Mapping, Postmap-QC, Variant Calling, Joint...
#   annotate
#   harmonize        Harmonize cromwell swarm script and input json...
#   joint_genotype   Joint Genotype cromwell swarm script and input json...
#   mapping          Mapping cromwell swarm script and input json preparation
#   mapping_report
#   postmap_qc       Post-mapping QC cromwell swarm script and input json...
#   premap_qc        Premap QC cromwell swarm script and input json...
#   premap_report    Premap report cromwell swarm script and input json...
#   variant_calling  Variant Calling cromwell swarm script and input json...

DEFAULT_SWARM_INPUT_TEMPLATE=['#SWARM --threads-per-process 8', '#SWARM --gb-per-process 25',  '#SWARM --time 4:00:00', '#SWARM --module cromwell,singularity', ]

DEFAULT_TEMPLATE_DIRECTORY = os.environ.get('GENCOMPASS_INPUT_TEMPLATE_DIRECTORY', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'templates'))
DEFAULT_PROJECT_PARAMETERS = os.environ.get('GENCOMPASS_PROJECT_PARAMETERS', os.path.join(DEFAULT_TEMPLATE_DIRECTORY, 'project_parameters.json'))
DEFAULT_PREMAP_INPUT_TEMPLATE =  os.path.join(DEFAULT_TEMPLATE_DIRECTORY, 'premap_qc.input_template.json')
DEFAULT_MAPPING_INPUT_TEMPLATE = os.path.join(DEFAULT_TEMPLATE_DIRECTORY, 'mapping_fq2bam.input_template.json')
DEFAULT_POSTMAP_QC_INPUT_TEMPLATE = os.path.join(DEFAULT_TEMPLATE_DIRECTORY, 'post_mapping_qc.input_template.json')
DEFAULT_VARIANT_CALLING_INPUT_TEMPLATE = os.path.join(DEFAULT_TEMPLATE_DIRECTORY, 'variant_calling.input_template.json')
DEFAULT_JOINT_GENOTYPE_INPUT_TEMPLATE = os.path.join(DEFAULT_TEMPLATE_DIRECTORY, 'joint_genotype.input_template.json')
DEFAULT_HARMONIZE_INPUT_TEMPLATE = os.path.join(DEFAULT_TEMPLATE_DIRECTORY, 'harmonize.input_template.json')
DEFAULT_PREMAP_REPORT_INPUT_TEMPLATE = os.path.join(DEFAULT_TEMPLATE_DIRECTORY, 'premap_report.input_template.json')
DEFAULT_MAPPING_REPORT_INPUT_TEMPLATE = os.path.join(DEFAULT_TEMPLATE_DIRECTORY, 'mapping_report.input_template.json')
DEFAULT_ANNOTATION_INPUT_TEMPLATE = os.path.join(DEFAULT_TEMPLATE_DIRECTORY, 'annotation.input_template.json')

@click.group(chain=True)
def cli():
    pass


def update_template_with_project_data(template: dict, project_dict: dict):
    for proj_key, proj_value in project_dict.items():
        for input_key, input_value in template.items():
            if '{{' + str(proj_key) + '}}' in str(input_value):
                input_value = input_value.replace('{{' + str(proj_key) + '}}', proj_value)
            template[input_key] = input_value
    return template

def get_sample_ids(project_parameters):
    if 'samples' in project_parameters.keys():
        return [s.strip() for s in open(project_parameters['samples'])]
    manifest = Manifest(project_parameters['manifest'])
    return manifest.get_sample_ids()



def build_cromwell_options(project_parameters, name):
    options = {
        "delete_intermediate_output_files": True,
        "final_workflow_outputs_dir": project_parameters['results_directory'],
        "use_relative_output_paths": True,
        "final_workflow_log_dir": project_parameters['logs_directory'],
        "final_call_logs_dir": project_parameters['logs_directory'],
    }
    options_dir = os.path.join(project_parameters['working_directory'], 'workflow_options')
    os.makedirs(options_dir, exist_ok=True)
    options_filename = os.path.join(options_dir, f"{name}_options.json")
    with open(options_filename, "w") as ofile:
        json.dump(options, ofile, indent=4)   


def create_fastq_files_table(fastq_files_path: str, manifest: Manifest, samples = []) :
    def get_sample_run_id(fastq_path: str) -> str:
        fastq_filename = os.path.basename(fastq_path)
        sample_id  = fastq_filename.split('_')[0]
        return sample_id
    fastq_files  = pd.read_csv(fastq_files_path, header=None, names=['Fastq Location'])
    fastq_files['Sample Run ID'] = fastq_files['Fastq Location'].apply(get_sample_run_id)

    fastq_files['Sample ID'] = fastq_files['Sample Run ID'].map(manifest.get_sample_run_id_map())
    #   = fastq_files.merge(manifest[['Sample ID', 'Sample Run ID']], how ='inner', on = 'Sample Run ID')
    fastq_files  = fastq_files [['Sample ID', 'Fastq Location']]
    if samples is not None and len(samples) > 0:
        fastq_files = fastq_files[fastq_files['Sample ID'].isin(samples)]
    return fastq_files 


def build_sample_level_inputs(project_parameters, input_template, name):
    input_params_dir= os.path.join(f"{project_parameters['working_directory']}","inputs", f"{name}_inputs")
    samples = get_sample_ids(project_parameters)
    os.makedirs(input_params_dir, exist_ok=True)
    for sample in samples:
        sample_input = input_template.copy()
        for key, value in sample_input.items():
            if '{{' + 'sample_id' + '}}' in str(value):
                sample_input[key] = value.replace('{{' + 'sample_id' + '}}', sample)
        with open(os.path.join(input_params_dir,f"{sample}_{name}_input.json"), "w") as ofile:
            json.dump(sample_input, ofile, indent=4)


def prepare_sample_level_files(project_parameters: dict, input_template: dict, name: str, wdl_filename):
    input_template = update_template_with_project_data(input_template, project_parameters)
    # copy wdl file to working directory
    wdl_abs_path = '/'.join((os.path.abspath(__file__)).split('/')[:-2]) + f'/workflows/{wdl_filename}'
    working_workflow_directory = os.path.join(project_parameters['working_directory'], 'workflows')
    os.makedirs(working_workflow_directory, exist_ok=True)

    shutil.copy(wdl_abs_path, working_workflow_directory)
    # create and save options file
    build_cromwell_options(project_parameters, name)
    # create and save sample level inputs
    build_sample_level_inputs(project_parameters, input_template, name)


def prepare_project_level_files(name, wdl_filename, project_parameters, input_template ):
    working_dir = project_parameters['working_directory'] if 'working_directory' in project_parameters.keys() else './'
    project_parameters['working_directory'] = working_dir
    build_cromwell_options(project_parameters, name)

    input_template = json.load(open(input_template))
    input_template = update_template_with_project_data(input_template, project_parameters)
    input_params_dir= os.path.join(f"{project_parameters['working_directory']}","inputs", f"{name}_inputs")
    os.makedirs(input_params_dir, exist_ok=True)

    # copy wdl file to working directory
    wdl_abs_path = '/'.join((os.path.abspath(__file__)).split('/')[:-2]) + f'/workflows/{wdl_filename}'
    working_workflow_directory = os.path.join(project_parameters['working_directory'], 'workflows')
    os.makedirs(working_workflow_directory, exist_ok=True)

    shutil.copy(wdl_abs_path, working_workflow_directory)
    input_json_path = os.path.join(input_params_dir, f'{name}.inputs.json')
    with open(input_json_path, 'w') as ofile:
            json.dump(input_template, ofile, indent=4)
    options_json_path = os.path.join(project_parameters['working_directory'], "workflow_options", f"{name}_options.json")

    cromwell_invocation = f'java -Dconfig.file={project_parameters["cromwell_config"]} -jar {project_parameters["cromwell_jar"]} run {os.path.join(working_workflow_directory, wdl_filename)} -o {options_json_path} -i {input_json_path}'
    singularity_cache = project_parameters['singularity_cachedir']
    swarm_log_dir = os.path.join(project_parameters["working_directory"], 'swarm_logs', name)
    swarm_parameters = DEFAULT_SWARM_INPUT_TEMPLATE.copy()
    swarm_parameters.append(f'#SWARM --logdir {swarm_log_dir}')
    swarm_parameters.append(f'#SWARM --sbatch "--export SINGULARITY_CACHEDIR={singularity_cache}"')
    with open(os.path.join(project_parameters['working_directory'], f"{project_parameters['project']}_{name}.swarm"), 'w') as ofile:
        ofile.write('\n'.join(swarm_parameters) + '\n')
        ofile.write(cromwell_invocation + '\n')
    

def create_sample_swarm_script(project_parameters, name, wdl_filename):
        input_params_dir= os.path.join(f"{project_parameters['working_directory']}","inputs", f"{name}_inputs")
        
        singularity_cache = project_parameters['singularity_cachedir']
        swarm_parameters = DEFAULT_SWARM_INPUT_TEMPLATE.copy()
        swarm_log_dir = os.path.join(project_parameters["working_directory"], 'swarm_logs', name)
        swarm_parameters.append(f'#SWARM --logdir {swarm_log_dir}')
        swarm_parameters.append(f'#SWARM --sbatch "--export SINGULARITY_CACHEDIR={singularity_cache}"')
        options_json_path = os.path.join(project_parameters['working_directory'], "workflow_options", f"{name}_options.json")
        working_workflow_path = os.path.join(project_parameters['working_directory'], 'workflows', wdl_filename)

        samples = get_sample_ids(project_parameters)

        with open(os.path.join(project_parameters['working_directory'], f"{project_parameters['project']}_{name}.swarm"), 'w') as ofile:
            ofile.write('\n'.join(swarm_parameters) + '\n')
            for sample in samples:
                sample = sample.strip()
                input_json = os.path.join(input_params_dir, f"{sample}_{name}_input.json")
                cromwell_invocation = f'java -Dconfig.file={project_parameters["cromwell_config"]} -jar {project_parameters["cromwell_jar"]} run {working_workflow_path} -o {options_json_path} -i {input_json}'
                ofile.write(cromwell_invocation + '\n')
        

@cli.command('premap_qc', help='Premap QC cromwell swarm script and input json preparation')
@click.option('--project_parameters', 
              help='Project Parameters - JSON file', 
              default=DEFAULT_PROJECT_PARAMETERS, 
              show_default=DEFAULT_PROJECT_PARAMETERS)
@click.option('--input_template', 
              help='Input Template - JSON file',
              default=DEFAULT_PREMAP_INPUT_TEMPLATE,
              show_default=DEFAULT_PREMAP_INPUT_TEMPLATE)
def premap(project_parameters,input_template):
    name='premap'
    wdl_filename = 'premap_qc.wdl'
    project_parameters = json.load(open(project_parameters))
    working_dir = project_parameters['working_directory'] if 'working_directory' in project_parameters.keys() else './'
    project_parameters['working_directory'] = working_dir
    input_template = json.load(open(input_template))
    prepare_sample_level_files(project_parameters, input_template, name, wdl_filename)
    create_sample_swarm_script(project_parameters, name, wdl_filename)


@cli.command('mapping', help='Mapping cromwell swarm script and input json preparation')
@click.option('--project_parameters', 
              help='Project Parameters - JSON file', 
              default=DEFAULT_PROJECT_PARAMETERS, 
              show_default=DEFAULT_PROJECT_PARAMETERS)
@click.option('--input_template', 
              help='Input Template - JSON file',
              default=DEFAULT_MAPPING_INPUT_TEMPLATE,
              show_default=DEFAULT_MAPPING_INPUT_TEMPLATE)
def mapping(project_parameters,input_template):

    wdl_filename = 'mapping.wdl'
    name='mapping'
    project_parameters = json.load(open(project_parameters))
    working_dir = project_parameters['working_directory'] if 'working_directory' in project_parameters.keys() else './'
    project_parameters['working_directory'] = working_dir
    input_template = json.load(open(input_template))

    manifest = Manifest(project_parameters['manifest'])
    fastq_table = create_fastq_files_table(project_parameters['fastq_files'], manifest) 

    def get_fastp_result_loc(sample: str, original_fastq_filename: str):
        lane_id = helpers.extract_sample_lane_id(original_fastq_filename)
        paired_end = helpers.extract_sample_paired_end(original_fastq_filename)
        fastp_fastq_path = os.path.join(project_parameters['results_directory'], f'premap_qc/fastp/{sample}/{lane_id}_{paired_end}_fastp.fastq.gz')
        return fastp_fastq_path
    fastq_table['fastp_fastq_path'] = fastq_table.apply(lambda row: get_fastp_result_loc(row['Sample ID'], row['Fastq Location']), axis = 1)

    input_template = update_template_with_project_data(input_template, project_parameters)
    # copy wdl file to working directory
    wdl_abs_path = '/'.join((os.path.abspath(__file__)).split('/')[:-2]) + f'/workflows/{wdl_filename}'
    working_workflow_directory = os.path.join(project_parameters['working_directory'], 'workflows')
    os.makedirs(working_workflow_directory, exist_ok=True)

    shutil.copy(wdl_abs_path, working_workflow_directory)
    # create and save options file
    build_cromwell_options(project_parameters, name)
    # create and save sample level inputs
    input_params_dir= os.path.join(f"{project_parameters['working_directory']}","inputs", f"{name}_inputs")
    os.makedirs(input_params_dir, exist_ok=True)
    samples = get_sample_ids(project_parameters)
    for sample in samples:
        sample = sample.strip()
        sample_input = input_template.copy()
        for key, value in sample_input.items():
            if '{{' + 'sample_id' + '}}' in str(value):
                sample_input[key] = value.replace('{{' + 'sample_id' + '}}', sample)

        fastp_fastq_files = fastq_table[fastq_table['Sample ID'] == sample]
        sample_fastq_files = list(fastp_fastq_files['fastp_fastq_path'].unique())
        sample_input['Mapping.sampleFastqFiles'] = sample_fastq_files

        with open(os.path.join(input_params_dir, f"{sample}_{name}_input.json"), "w") as ofile:
            json.dump(sample_input, ofile, indent=4)
    create_sample_swarm_script(project_parameters, name, wdl_filename)


@cli.command('postmap_qc', help='Post-mapping QC cromwell swarm script and input json preparation')
@click.option('--project_parameters', 
              help='Project Parameters - JSON file', 
              default=DEFAULT_PROJECT_PARAMETERS, 
              show_default=DEFAULT_PROJECT_PARAMETERS)
@click.option('--input_template', 
              help='Input Template - JSON file',
              default=DEFAULT_POSTMAP_QC_INPUT_TEMPLATE,
              show_default=DEFAULT_POSTMAP_QC_INPUT_TEMPLATE)
def postmap_qc(project_parameters,input_template):
    name = 'postmap_qc'
    wdl_filename='mapping.wdl'
    project_parameters = json.load(open(project_parameters))
    working_dir = project_parameters['working_directory'] if 'working_directory' in project_parameters.keys() else './'
    project_parameters['working_directory'] = working_dir
    input_template = json.load(open(input_template))
    prepare_sample_level_files(project_parameters, input_template, name, wdl_filename)
    create_sample_swarm_script(project_parameters, name, wdl_filename)


@cli.command('variant_calling', help='Variant Calling cromwell swarm script and input json preparation')
@click.option('--project_parameters', 
              help='Project Parameters - JSON file', 
              default=DEFAULT_PROJECT_PARAMETERS, 
              show_default=DEFAULT_PROJECT_PARAMETERS)
@click.option('--input_template', 
              help='Input Template - JSON file',
              default=DEFAULT_VARIANT_CALLING_INPUT_TEMPLATE,
              show_default=DEFAULT_VARIANT_CALLING_INPUT_TEMPLATE)
def variant_calling(project_parameters,input_template):
    name = 'variant_calling'
    wdl_filename='variant_calling.wdl'
    project_parameters = json.load(open(project_parameters))
    working_dir = project_parameters['working_directory'] if 'working_directory' in project_parameters.keys() else './'
    project_parameters['working_directory'] = working_dir
    input_template = json.load(open(input_template))
    prepare_sample_level_files(project_parameters, input_template, name, wdl_filename)
    create_sample_swarm_script(project_parameters, name, wdl_filename)

def get_joint_genotype_runtime_request(sample_size, multiplier=3):
    runtime_minutes = sample_size * multiplier + 60
    runtime_hours = runtime_minutes // 60
    runtime_minutes_remain = runtime_minutes % 60
    runtime_days = runtime_minutes // (60 * 24)
    runtime_hours = (runtime_minutes % (60 * 24)) // 60
    runtime_request = f'{runtime_minutes_remain:2}:00'
    runtime_request = f'{runtime_hours}:{runtime_request}' if runtime_hours > 0 else runtime_request
    runtime_request = f'{runtime_days}-{runtime_request}' if runtime_days > 0 else runtime_request
    return runtime_request


@cli.command('joint_genotype', help='Joint Genotype cromwell swarm script and input json preparation')
@click.option('--project_parameters', 
              help='Project Parameters - JSON file', 
              default=DEFAULT_PROJECT_PARAMETERS, 
              show_default=DEFAULT_PROJECT_PARAMETERS)
@click.option('--input_template', 
              help='Input Template - JSON file',
              default=DEFAULT_JOINT_GENOTYPE_INPUT_TEMPLATE,
              show_default=DEFAULT_JOINT_GENOTYPE_INPUT_TEMPLATE)
def joint_genptype(project_parameters,input_template):
    wdl_filename='joint_genotype.wdl'
    name='joint_genotype'
    project_parameters = json.load(open(project_parameters))
    working_dir = project_parameters['working_directory'] if 'working_directory' in project_parameters.keys() else './'
    project_parameters['working_directory'] = working_dir
    input_template = json.load(open(input_template))

    samples = get_sample_ids(project_parameters)

    input_template = update_template_with_project_data(input_template, project_parameters)
    input_template["JointGenotype.glnexusRuntimeAttributes"]["runtimeMinutes"] = len(samples) * 3 + 60
    # copy wdl file to working directory
    wdl_abs_path = '/'.join((os.path.abspath(__file__)).split('/')[:-2]) + f'/workflows/{wdl_filename}'
    working_workflow_directory = os.path.join(project_parameters['working_directory'], 'workflows')
    os.makedirs(working_workflow_directory, exist_ok=True)

    shutil.copy(wdl_abs_path, working_workflow_directory)
    input_params_dir= os.path.join(f"{project_parameters['working_directory']}","inputs", f"{name}_inputs")
    
    os.makedirs(input_params_dir, exist_ok=True)
    shutil.copy('/'.join((os.path.abspath(__file__)).split('/')[:-2]) + f'/config/strelka2_glnexus.yml', input_params_dir)

    glnexus_config_map = {'haplotypecaller': 'gatk',
                              'deepvariant': 'DeepVariant',
                              'strelka2': os.path.join(input_params_dir, 'strelka2_glnexus.yml') }
    # create and save options file
    build_cromwell_options(project_parameters, name)
    # create and save sample level inputs
    for caller, glnexus_config in glnexus_config_map.items():
        variant_fof = os.path.join(input_params_dir ,f'{caller}_vcf_files.txt')
        caller_input = input_template.copy()
        for input_key, input_value in caller_input.items():
            if '{{caller}}' in str(input_value):
                caller_input[input_key] = input_value.replace('{{caller}}', caller)
            if '{{glnexus_config}}' in str(input_value):
                caller_input[input_key] = input_value.replace('{{glnexus_config}}', glnexus_config)
            if '{{variants_file}}' in str(input_value):
                caller_input[input_key] = input_value.replace('{{variants_file}}', variant_fof)
        
        
        
        with open(os.path.join(input_params_dir, f'joint_genotype.{caller}_input_parameters.json'), 'w') as ofile:
            json.dump(caller_input, ofile, indent=4)
        
        output_extenstion_map = {'deepvariant' : '.deepvariant.g.vcf.gz', 
                                 'haplotypecaller' : '.haplotypecaller.g.vcf.gz',
                                 'strelka2': '.strelka.genome.vcf.gz'}
        with open(variant_fof, 'w') as ofile:
            for sample in samples:
                ofile.write(os.path.join(project_parameters['results_directory'], caller, sample, f'{sample}{output_extenstion_map[caller]}\n'))
            
    
    
    singularity_cache = project_parameters['singularity_cachedir']
    runtime_request = get_joint_genotype_runtime_request(len(samples), 3)
    swarm_parameters = ['#SWARM --threads-per-process 8', '#SWARM --gb-per-process 25',  f'#SWARM --time {runtime_request}', '#SWARM --module cromwell,singularity', ]

    swarm_log_dir = os.path.join(project_parameters["working_directory"], 'swarm_logs', name)
    swarm_parameters.append(f'#SWARM --logdir {swarm_log_dir}')
    swarm_parameters.append(f'#SWARM --sbatch "--export SINGULARITY_CACHEDIR={singularity_cache}"')
    
    options_json_path = os.path.join(project_parameters['working_directory'], "workflow_options", f"{name}_options.json")
    workflow_path=os.path.join(project_parameters["working_directory"],"workflows", wdl_filename)
    with open(os.path.join(project_parameters['working_directory'], f"{project_parameters['project']}_{name}.swarm"), 'w') as ofile:
        ofile.write('\n'.join(swarm_parameters) + '\n')
        for caller in glnexus_config_map.keys():
            input_json_path = os.path.join(input_params_dir, f'joint_genotype.{caller}_input_parameters.json')
            cromwell_invocation = f'java -Dconfig.file={project_parameters["cromwell_config"]} -jar {project_parameters["cromwell_jar"]} run {workflow_path} -o {options_json_path} -i {input_json_path}'
            ofile.write(cromwell_invocation + '\n')


@cli.command('harmonize', help='Harmonize cromwell swarm script and input json preparation')
@click.option('--project_parameters', 
              help='Project Parameters - JSON file', 
              default=DEFAULT_PROJECT_PARAMETERS, 
              show_default=DEFAULT_PROJECT_PARAMETERS)
@click.option('--joint_genotype_input_template', 
              help='Joint genotype input template - JSON file',
              default=DEFAULT_JOINT_GENOTYPE_INPUT_TEMPLATE,
              show_default=DEFAULT_JOINT_GENOTYPE_INPUT_TEMPLATE
              )
@click.option('--harmonize_input_template', 
              help='Harmonize input template - JSON file',
              default=DEFAULT_HARMONIZE_INPUT_TEMPLATE,
              show_default=DEFAULT_HARMONIZE_INPUT_TEMPLATE
              )
def harmonize(project_parameters,joint_genotype_input_template, harmonize_input_template):
    name='harmonize'
    wdl_filename='harmonize.wdl'
    project_parameters = json.load(open(project_parameters))
    working_dir = project_parameters['working_directory'] if 'working_directory' in project_parameters.keys() else './'
    project_parameters['working_directory'] = working_dir
    build_cromwell_options(project_parameters, name)

    harmonize_input_template = json.load(open(harmonize_input_template))
    harmonize_input_template = update_template_with_project_data(harmonize_input_template, project_parameters)
    joint_genotype_input_template = json.load(open(joint_genotype_input_template))
    joint_genotype_input_template = update_template_with_project_data(joint_genotype_input_template, project_parameters)
    interval_names = get_interval_names(bedfile=joint_genotype_input_template['JointGenotype.intervalBED'], nucleotides_per_bin=joint_genotype_input_template['JointGenotype.nucleotidesPerBin'])
    
    input_params_dir= os.path.join(f"{project_parameters['working_directory']}", "inputs",f"{name}_inputs")
    os.makedirs(input_params_dir, exist_ok=True)

    # copy wdl file to working directory
    wdl_abs_path = '/'.join((os.path.abspath(__file__)).split('/')[:-2]) + f'/workflows/{wdl_filename}'
    working_workflow_directory = os.path.join(project_parameters['working_directory'], 'workflows')
    os.makedirs(working_workflow_directory, exist_ok=True)

    shutil.copy(wdl_abs_path, working_workflow_directory)
    
    for bin in interval_names:
        bin_harmonize_input = harmonize_input_template.copy()
        for input_key, input_value in bin_harmonize_input.items():
            if '{{bin}}' in str(input_value):
                bin_harmonize_input[input_key] = input_value.replace('{{bin}}', bin)
        print(bin)
        
        with open(os.path.join(input_params_dir, f'harmonize.{bin}_inputs.json'), 'w') as ofile:
            json.dump(bin_harmonize_input, ofile, indent=4)

    singularity_cache = project_parameters['singularity_cachedir']
    swarm_parameters = DEFAULT_SWARM_INPUT_TEMPLATE.copy()
    swarm_log_dir = os.path.join(project_parameters["working_directory"], 'swarm_logs', name)
    swarm_parameters.append(f'#SWARM --logdir {swarm_log_dir}')
    swarm_parameters.append(f'#SWARM --sbatch "--export SINGULARITY_CACHEDIR={singularity_cache}"')
    
    options_json_path = os.path.join(project_parameters['working_directory'], "workflow_options", f"{name}_options.json")
    workflow_path = os.path.join(project_parameters["working_directory"],"workflows", wdl_filename)
    with open(os.path.join(project_parameters['working_directory'], f"{project_parameters['project']}_{name}.swarm"), 'w') as ofile:
        ofile.write('\n'.join(swarm_parameters) + '\n')
        for bin in interval_names:
            input_json_path = os.path.join(input_params_dir, f'harmonize.{bin}_inputs.json')
            cromwell_invocation = f'java -Dconfig.file={project_parameters["cromwell_config"]} -jar {project_parameters["cromwell_jar"]} run {workflow_path} -o {options_json_path} -i {input_json_path}'
            ofile.write(cromwell_invocation + '\n')

@cli.command('premap_report', help='Premap report cromwell swarm script and input json preparation')
@click.option('--project_parameters', 
              help='Project Parameters - JSON file', 
              default=DEFAULT_PROJECT_PARAMETERS, 
              show_default=DEFAULT_PROJECT_PARAMETERS)
@click.option('--input_template', 
              help='Input Template - JSON file',
              default=DEFAULT_PREMAP_REPORT_INPUT_TEMPLATE,
              show_default=DEFAULT_PREMAP_REPORT_INPUT_TEMPLATE)
def premap_report(project_parameters,input_template):
    name='premap_report'
    wdl_filename='premap_qc_report.wdl'
    project_parameters = json.load(open(project_parameters))
    prepare_project_level_files(name, wdl_filename, project_parameters, input_template )
    
@cli.command('mapping_report')
@click.option('--project_parameters', 
              help='Project Parameters - JSON file', 
              default=DEFAULT_PROJECT_PARAMETERS, 
              show_default=DEFAULT_PROJECT_PARAMETERS)
@click.option('--input_template', 
              help='Input Template - JSON file',
              default=DEFAULT_MAPPING_REPORT_INPUT_TEMPLATE,
              show_default=DEFAULT_MAPPING_REPORT_INPUT_TEMPLATE)
def mapping_report(project_parameters,input_template):
    name='mapping_report'
    wdl_filename='mapping_report.wdl'
    project_parameters = json.load(open(project_parameters))
    prepare_project_level_files(name, wdl_filename, project_parameters, input_template )

@cli.command('annotate')
@click.option('--project_parameters', 
              help='Project Parameters - JSON file', 
              default=DEFAULT_PROJECT_PARAMETERS, 
              show_default=DEFAULT_PROJECT_PARAMETERS)
@click.option('--input_template', 
              help='Input Template - JSON file',
              default=DEFAULT_ANNOTATION_INPUT_TEMPLATE,
              show_default=DEFAULT_ANNOTATION_INPUT_TEMPLATE)
def annotate(project_parameters,input_template):
    pass

@cli.command('all', help='Premap, Mapping, Postmap-QC, Variant Calling, Joint Genotype, Harmonize, Annotate, Premap-QC Report, Mapping Report')
@click.option('--project_parameters', 
              help='Project Parameters - JSON file', 
              default=DEFAULT_PROJECT_PARAMETERS, 
              show_default=DEFAULT_PROJECT_PARAMETERS)
@click.option('--input_template_directory', 
              help='Input template directory',
              default=DEFAULT_TEMPLATE_DIRECTORY,
              show_default=DEFAULT_TEMPLATE_DIRECTORY)
@click.pass_context
def all(ctx,project_parameters, input_template_directory):
    premap_input_template =  os.path.join(input_template_directory, 'premap_qc.input_template.json')
    mapping_fq2bam_input_template = os.path.join(input_template_directory, 'mapping_fq2bam.input_template.json')
    postmap_qc_input_template = os.path.join(input_template_directory, 'post_mapping_qc.input_template.json')
    variant_calling_input_template = os.path.join(input_template_directory, 'variant_calling.input_template.json')
    joint_genptype_input_template = os.path.join(input_template_directory, 'joint_genotype.input_template.json')
    harmonize_input_template = os.path.join(input_template_directory, 'harmonize.input_template.json')
    premap_report_input_template = os.path.join(input_template_directory, 'premap_report.input_template.json')
    mapping_report_input_template = os.path.join(input_template_directory, 'mapping_report.input_template.json')
    annotation_input_template = os.path.join(input_template_directory, 'annotation.input_template.json')
    ctx.invoke(premap,project_parameters=project_parameters, input_template=premap_input_template)
    ctx.invoke(mapping,project_parameters=project_parameters, input_template=mapping_fq2bam_input_template)

    ctx.invoke(postmap_qc,project_parameters=project_parameters, input_template=postmap_qc_input_template)
    ctx.invoke(variant_calling,project_parameters=project_parameters, input_template=variant_calling_input_template)
    ctx.invoke(joint_genptype,project_parameters=project_parameters, input_template=joint_genptype_input_template)
    ctx.invoke(harmonize,project_parameters=project_parameters, joint_genotype_input_template=joint_genptype_input_template, harmonize_input_template=harmonize_input_template)
    ctx.invoke(premap_report,project_parameters=project_parameters, input_template=premap_report_input_template)
    ctx.invoke(mapping_report,project_parameters=project_parameters, input_template=mapping_report_input_template)
    ctx.invoke(annotate,project_parameters=project_parameters, input_template=annotation_input_template)



        


if __name__ == '__main__':
    cli()



