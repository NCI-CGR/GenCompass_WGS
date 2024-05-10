WORKFLOW_RESULTS_DIR=$1
GENCOMPASS=$2
PROJECT=$3
OUTPUT_DIR=$4


mkdir -p $OUTPUT_DIR
mkdir $OUTPUT_DIR/workflow_files
#-------------------------------------
# Sequencing Artifact Reports
#   - pre
#   - bait
#-------------------------------------

echo "Generating Sequencing Artifact Reports"
ls $WORKFLOW_RESULTS_DIR/mapping_qc/collectmultiplemetrics/*/sequencingArtifact.pre_adapter_summary_metrics.txt > $OUTPUT_DIR/workflow_files/sequencingArtifact.pre_adapter_summary_metrics.FOF.txt
ls $WORKFLOW_RESULTS_DIR/mapping_qc/collectmultiplemetrics/*/sequencingArtifact.bait_bias_summary_metrics.txt > $OUTPUT_DIR/workflow_files/sequencingArtifact.bait_bias_summary_metrics.FOF.txt

python $GENCOMPASS/scripts/sequencing_artifact_reports.py \
    --seqart_pre_list $OUTPUT_DIR/workflow_files/sequencingArtifact.pre_adapter_summary_metrics.FOF.txt \
    --seqart_bait_list $OUTPUT_DIR/workflow_files/sequencingArtifact.bait_bias_summary_metrics.FOF.txt \
    --project $PROJECT \
    -odir $OUTPUT_DIR


#-------------------------------------
# BamMetrics Report
#-------------------------------------
echo "Generating Bammetrics Report"
ls $WORKFLOW_RESULTS_DIR/mapping_qc/bammetrics/*.bammetrics.txt > $OUTPUT_DIR/workflow_files/bammetrics.FOF.txt

python $GENCOMPASS/scripts/bammetrics_report.py \
    --bammetrics_list $OUTPUT_DIR/workflow_files/bammetrics.FOF.txt \
    -odir $OUTPUT_DIR \
    --project $PROJECT


#-------------------------------------
# MultiQC report
#-------------------------------------
echo "Running MultiQC"
multiqc \
    -o $OUTPUT_DIR \
    -n mapping.multiqc \
    -t Mapping-MultiQC-Report \
    ${WORKFLOW_RESULTS_DIR}/mapping_qc/collectmultiplemetrics \
    ${WORKFLOW_RESULTS_DIR}/mapping_qc/verifybamid \
    ${WORKFLOW_RESULTS_DIR}/mapping_qc/somalier

#-------------------------------------
# Post Mapping report
#-------------------------------------
# echo "Running MultiQC on collectmultiplemetrics"
# multiqc \
#     -o $OUTPUT_DIR/multiqc \
#     -n collectmultiple_metrics.multiqc \
#     ${WORKFLOW_RESULTS_DIR}/mapping_qc/collectmultiplemetrics

MULTIQC_JSON=$OUTPUT_DIR/mapping.multiqc_data/multiqc_data.json

echo "Generating Post Mapping Report"
python $GENCOMPASS/scripts/post_mapping_report.py \
    --cmm_multiqc_json $MULTIQC_JSON \
    --bammetrics_report $OUTPUT_DIR/$PROJECT.bammetrics_report.tsv \
    --project $PROJECT \
    -odir $OUTPUT_DIR


#-------------------------------------
# VerifyBamID2 Report
#-------------------------------------
echo "Generating VerifyBamID2 Report"
ls $WORKFLOW_RESULTS_DIR/mapping_qc/verifybamid/*/*.selfSM > $OUTPUT_DIR/workflow_files/verifybamid.selfSM.FOF.txt
python $GENCOMPASS/scripts/verifybamid2_report.py \
    --selfsm_list $OUTPUT_DIR/workflow_files/verifybamid.selfSM.FOF.txt \
    -odir $OUTPUT_DIR \
    --project $PROJECT

#-------------------------------------
# Somalier Sex and Relatedness Reports
#-------------------------------------
echo "Generating Somalier Sex and Relatedness Reports"
python $GENCOMPASS/scripts/somalier_report.py \
    --somalier_samples $WORKFLOW_RESULTS_DIR/mapping_qc/somalier/relate/$PROJECT.samples.tsv \
    --somalier_pairs $WORKFLOW_RESULTS_DIR/mapping_qc/somalier/relate/$PROJECT.pairs.tsv \
    -odir $OUTPUT_DIR \
    --project $PROJECT \
    --ibs2_cutoff 5000 \
    --relatedness_cutoff 0.0

#-------------------------------------
# Kraken2 Report
#-------------------------------------
echo "Generating Kraken2 Report"
python $GENCOMPASS/scripts/kraken_report.py \
    --kraken_reports ${WORKFLOW_RESULTS_DIR}/mapping_qc/kraken2/*/*.kraken.mpaStyle.report \
    --kraken_outputs ${WORKFLOW_RESULTS_DIR}/mapping_qc/kraken2/*/*.kraken.output \
    --project $PROJECT \
    -odir $OUTPUT_DIR