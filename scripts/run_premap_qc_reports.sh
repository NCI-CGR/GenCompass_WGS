WORKFLOW_RESULTS_DIR=$1
MANIFEST=$2
GENCOMPASS=$3
PROJECT=$4
SEQUENCING_TYPE=$5
OUTPUT_DIR=$6

mkdir -p $OUTPUT_DIR

# fastp report
echo "Generating fastp report"
ls $WORKFLOW_RESULTS_DIR/premap_qc/fastp/*/*.json > $OUTPUT_DIR/fastp_json_files.txt
python $GENCOMPASS/scripts/fastp_report.py \
--fastp_json_list $OUTPUT_DIR/fastp_json_files.txt \
--manifest $MANIFEST \
--project_name $PROJECT \
--sequencing_type $SEQUENCING_TYPE \
--output_directory $OUTPUT_DIR

# fastq_screen report
echo "Generating fastq screen report"
ls $WORKFLOW_RESULTS_DIR/premap_qc/fastq_screen/*/*.txt > $OUTPUT_DIR/fastq_screen_files.txt
python $GENCOMPASS/scripts/fastq_screen_report.py \
--fastq_screen_list $OUTPUT_DIR/fastq_screen_files.txt \
--project_name $PROJECT \
--output_directory $OUTPUT_DIR

echo "Generating PremapQC multiqc report"
multiqc \
    --title PremapQC \
    --filename PremapQC_multiqc_report.html \
    --outdir $OUTPUT_DIR \
    ${WORKFLOW_RESULTS_DIR}/premap_qc \








