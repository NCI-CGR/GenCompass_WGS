# GenCompass Batch Builder - running at scale with cromwell and swarm

The `prepare_gencompass.py` script is used to build input json files and swarm scripts 
for the various workflows included in GenCompass. The script has several commands that can be chained together

## Before using

1. Create a manifest file
2. Create a sample list file
3. Copy `batch_builder/templates/project_parameters.json` to working directory, update parameters with project specific data
4. (Optional) Update input templates with project specific parameters if different from default
   1. Update docker if not using cgrlab DockerHub images
   2. Update task resource requirements if needed
5. (Optional) Set environmental variable `GENCOMPASS_PROJECT_PARAMETERS` (this can be also set with a command line argument)

### Example 1. Generating all workflow files

```bash
# Provide project parameters json in CLI. Optinally provide an input template directory if needed.
python GenCompass_WGS/batch_builder/prepare_gencompass.py \
all \
--project_parameters ./my_project_parameters.json

```

### Example 2. Generate single workflow inputs
```bash
python Gencompass/batch_builder/prepare_gencompass.py \
premap_qc \
--project_parameters ./my_project_parameters.json \
--input_template ./project_input_templates/premap_qc.input_template.json
```
Single workflows available to generate: `premap_qc`, `premap_report`, `mapping`, `postmap_qc`, `mapping_report`, `variant_calling`, `joint_genotype`, `harmonize`


## Input Files

### Manifest
- The manifest file can be a csv or xlsx file. 
- Required columns: `Sample ID`, `Flowcell`
- Optional column : `Sample Run ID`

Sample ID|Flowcell|Sample Run ID
---------|--------|--------------
EXAMPLE-SAMPLE-136|B00FLOWCELL|SAMPLE_RUN-ID-6292
EXAMPLE-SAMPLE-003|B00FLOWCELL|SAMPLE_RUN-ID-5012
EXAMPLE-SAMPLE-082|A00FLOWCELL|SAMPLE_RUN-ID-4014
EXAMPLE-SAMPLE-145|A00FLOWCELL|SAMPLE_RUN-ID-4751
EXAMPLE-SAMPLE-159|A00FLOWCELL|SAMPLE_RUN-ID-3034
EXAMPLE-SAMPLE-182|A00FLOWCELL|SAMPLE_RUN-ID-7312

### Sample List (optional)
**example_samples.txt**
```
EXAMPLE-SAMPLE-136
EXAMPLE-SAMPLE-003
```

### Fastq Files
**example_fastq_files.txt**
```
/path/to/fastq/SAMPLE_RUN-ID-6292_S2_L001_R1_001.fastq.gz
/path/to/fastq/SAMPLE_RUN-ID-6292_S2_L001_R2_001.fastq.gz
/path/to/fastq/SAMPLE_RUN-ID-6292_S2_L002_R1_001.fastq.gz
/path/to/fastq/SAMPLE_RUN-ID-6292_S2_L002_R2_001.fastq.gz
/path/to/fastq/SAMPLE_RUN-ID-5012_S2_L001_R1_001.fastq.gz
/path/to/fastq/SAMPLE_RUN-ID-5012_S2_L001_R2_001.fastq.gz
/path/to/fastq/SAMPLE_RUN-ID-5012_S2_L002_R1_001.fastq.gz
/path/to/fastq/SAMPLE_RUN-ID-5012_S2_L002_R2_001.fastq.gz
```

### Project Parameters JSON

> Note: If `samples` is not provided then all samples in manifest will be used.

**project_parameters.json**
```json
{
    "project": "ExampleProject",
    "manifest": "Example_Manifest.csv",
    "fastq_files" : "example_fastq_files.txt",
    "reference_bundle": "/path/to/ReferenceBundle",
    "samples" : "example_samples.txt",
    "results_directory": "workflow_results",
    "logs_directory" :"workflow_logs",
    "working_directory" : "./",
    "cromwell_jar": "$CROMWELL_JAR",
    "cromwell_config": "$CROMWELL_CONFIG",
    "singularity_cachedir": "./singularity_cache"
}
```

## Running the workflows

The prepare_gencompass.py script generates a swarm file and all required input files required to run each workflow. Each swarm file will start with the `project` parameter defined in the `project_parameters.json` file. Each workflow is run in order 

### Data Generation Workflows
1. ExampleProject_premap.swarm
2. ExampleProject_mapping.swarm
3. ExampleProject_postmap_qc.swarm
4. ExampleProject_variant_calling.swarm
5. ExampleProject_joint_genotype.swarm
6. ExampleProject_harmonize.swarm

### Report Workflows
1. ExampleProject_premap_report.swarm
2. ExampleProject_mapping_report.swarm
