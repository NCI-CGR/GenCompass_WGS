## GenCompass: Germline Ensemble Calling of Mutations with Parabricks Accelerated Software Suite

---

### AUTHORS

Ben Jordan and Laura Egolf

Technical Lead: Komal Jain

---

## Table of Contents

1. [Running GenCompass](#running-gencompass)
2. [Input Files](#input-files)

## Running GenCompass
The `batch_builder/prepare_gencompass.py` script is used to build input json files and swarm scripts 
for the various workflows included in GenCompass. The script has several commands that can be chained together

### Before using

1. Create a [manifest file](#manifest)
2. Create a [list of fastq files](#fastq-files)
3. (Optional) Create a [sample list file](#sample-list-optional)
4. Copy `batch_builder/templates/project_parameters.json` to working directory, update parameters with project specific data
5. (Optional) Update input templates with project specific parameters if different from default
   1. Update docker if not using cgrlab DockerHub images
   2. Update task resource requirements if needed
6. (Optional) Set environmental variable `GENCOMPASS_PROJECT_PARAMETERS` (this can be also set with a command line argument)

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

`example_manifest.csv`
Sample ID|Flowcell|Sample Run ID
---------|--------|--------------
EXAMPLE-SAMPLE-136|B00FLOWCELL|SAMPLE_RUN-ID-6292
EXAMPLE-SAMPLE-003|B00FLOWCELL|SAMPLE_RUN-ID-5012
EXAMPLE-SAMPLE-082|A00FLOWCELL|SAMPLE_RUN-ID-4014
EXAMPLE-SAMPLE-082|A00FLOWCELL|SAMPLE_RUN-ID-4751


### Sample List (optional)
`example_samples.txt`
```
EXAMPLE-SAMPLE-136
EXAMPLE-SAMPLE-003
```

### Fastq Files
`example_fastq_files.txt`
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
