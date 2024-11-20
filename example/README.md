# GenCompass WGS Example

The following example has been created in order to demonstrate how to set up and run GenCompass.

## Set up

1. Download and copy reference bundle to example directory

## Generate scripts

```bash
python ../batch_builder/prepare_gencompass.py all --project_parameters example_project_parameters.json
```

The batch builder generates swarm files, used on slurm systems to submit many jobs at once. Each line is a command which can be run individually.

| File | Order | Description |
|------|-------|-------------|
example_premap.swarm | 1 | Premap QC |
example_premap_report.swarm | 2 | Premap QC report |
example_mapping.swarm | 3 | Mapping with fq2bam |
example_postmap_qc.swarm | 4 | Postmap QC |
example_mapping_report.swarm | 5 | Mapping Report |
example_variant_calling.swarm | 6 | Variant Calling |
example_joint_genotype.swarm | 7 | Joint Genotype |
example_harmonize.swarm | 8 | Harmonize |




**inputs folder structure**

```bash
inputs
├── harmonize_inputs
│   ├── harmonize.chr1.0_inputs.json
│   ├── harmonize.chr1.1_inputs.json
│   ├── harmonize.chr10.0_inputs.json
│   ├── harmonize.chr10.1_inputs.json
│   ├── harmonize.chr11.0_inputs.json
│   ├── harmonize.chr11.1_inputs.json
│   ├── harmonize.chr12_inputs.json
│   ├── harmonize.chr13_inputs.json
│   ├── harmonize.chr14_inputs.json
│   ├── harmonize.chr15_inputs.json
│   ├── harmonize.chr16_inputs.json
│   ├── harmonize.chr17_inputs.json
│   ├── harmonize.chr18_inputs.json
│   ├── harmonize.chr19_inputs.json
│   ├── harmonize.chr2.0_inputs.json
│   ├── harmonize.chr2.1_inputs.json
│   ├── harmonize.chr20_inputs.json
│   ├── harmonize.chr21_inputs.json
│   ├── harmonize.chr22_inputs.json
│   ├── harmonize.chr3.0_inputs.json
│   ├── harmonize.chr3.1_inputs.json
│   ├── harmonize.chr4.0_inputs.json
│   ├── harmonize.chr4.1_inputs.json
│   ├── harmonize.chr5.0_inputs.json
│   ├── harmonize.chr5.1_inputs.json
│   ├── harmonize.chr6.0_inputs.json
│   ├── harmonize.chr6.1_inputs.json
│   ├── harmonize.chr7.0_inputs.json
│   ├── harmonize.chr7.1_inputs.json
│   ├── harmonize.chr8.0_inputs.json
│   ├── harmonize.chr8.1_inputs.json
│   ├── harmonize.chr9.0_inputs.json
│   ├── harmonize.chr9.1_inputs.json
│   ├── harmonize.chrX.0_inputs.json
│   ├── harmonize.chrX.1_inputs.json
│   └── harmonize.chrY_inputs.json
├── joint_genotype_inputs
│   ├── deepvariant_vcf_files.txt
│   ├── haplotypecaller_vcf_files.txt
│   ├── joint_genotype.deepvariant_input_parameters.json
│   ├── joint_genotype.haplotypecaller_input_parameters.json
│   ├── joint_genotype.strelka2_input_parameters.json
│   ├── strelka2_glnexus.yml
│   └── strelka2_vcf_files.txt
├── mapping_inputs
│   ├── SM123456_mapping_input.json
│   └── SM246810_mapping_input.json
├── mapping_report_inputs
│   └── mapping_report.inputs.json
├── postmap_qc_inputs
│   ├── SM123456_postmap_qc_input.json
│   └── SM246810_postmap_qc_input.json
├── premap_inputs
│   ├── SM123456_premap_input.json
│   └── SM246810_premap_input.json
├── premap_report_inputs
│   └── premap_report.inputs.json
└── variant_calling_inputs
    ├── SM123456_variant_calling_input.json
    └── SM246810_variant_calling_input.json
```