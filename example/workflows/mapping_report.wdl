version 1.0

workflow MappingReport {
    input {

        File? sampleList
        String? workflowResultsDirectory

        File? somalierExtractsFOF
        File? sequencingArtifactPreAdapterSummaryMetricsFOF
        File? sequencingArtifactBaitBiasSummaryMetricsFOF
        File? bammetricsFOF 
        File? collectmultipleMetricsFOF # collectmultiplemetrics .txt files other than sequencing artifact files
        File? verifybamidFOF # verifybamid self.SM files
        File? krakenReportFOF # kraken 
        File? krakenOutputFOF

        String project
        String gencompassDocker
        String somalierDocker
        String multiqcDocker

        Boolean runSomalierReport=true
        Boolean runSequencingArtifactReport=true
        Boolean runMultiQCReport=true
        Boolean runParabricksReport=true
        Boolean runVerifyBamIdReport=true
        Boolean runKrakenReport=true
    }
    

    # call extractDate
    
    if( defined(sampleList) && defined(workflowResultsDirectory)) {
        call buildFileList {
            input:
                sampleList = select_first([sampleList]),
                workflowResultsDirectory = select_first([workflowResultsDirectory]),
                docker=gencompassDocker
        }
    }

    Array[File] verifyBamIdSelfSM = read_lines(select_first([verifybamidFOF, buildFileList.verifybamidSelfSmFiles]))

    if(runSomalierReport) {
        call somalierRelate {
            input:
                extracted = read_lines(select_first([somalierExtractsFOF, buildFileList.somalierExtractFiles])),
                project = project,
                docker = somalierDocker
        }
        call somalierSexAndRelatednessReport {
            input:
                somalierPairs = somalierRelate.relatednessPairs,
                somalierSamples = somalierRelate.samplesTSV,
                docker = gencompassDocker,
                project = project,
                # date = extractDate.date,
                oDir = "mapping_qc/reports"
            
        }
    }
        
    if(runSequencingArtifactReport){
        call sequencingArtifactReport {
            input:
                preAdapterSummaryMetrics = read_lines(select_first([sequencingArtifactPreAdapterSummaryMetricsFOF, buildFileList.cmmPreAdapterFiles])),
                baitBiasSummaryMetrics = read_lines(select_first([sequencingArtifactBaitBiasSummaryMetricsFOF, buildFileList.cmmBaitBiasFiles])),
                docker = gencompassDocker,
                project = project,
                # date = extractDate.date,
                oDir = "mapping_qc/reports"
        }

    }
    

    if(runParabricksReport){
        call bammetricsReport {
            input:
                bammetricsIndividualReports = read_lines(select_first([bammetricsFOF, buildFileList.bammetricsFiles])),
                docker = gencompassDocker,
                project = project,
                oDir = "mapping_qc/reports"
        }
    }
    
    if(runMultiQCReport){
        call multiqcReport {
            input:
                inputData = flatten([read_lines(select_first([collectmultipleMetricsFOF, buildFileList.cmmSummaryFiles])), verifyBamIdSelfSM, select_all([somalierRelate.relatednessPairs, somalierRelate.samplesTSV])]),
                docker = multiqcDocker,
                title = "~{project}-Mapping-MultiQC-Report",
                filename = "~{project}_mapping_multiqc",
                oDir = "mapping_qc/reports"
        }
    }
    
    if(runParabricksReport && runMultiQCReport){
        call mappingParabricksReport {
        input:
            bammetricsReportTable = select_first([bammetricsReport.bammetricsReportTable]),
            multiqcJSON = select_first([multiqcReport.multiqcDataJson]),
            docker = gencompassDocker,
            project = project,
            # date = extractDate.date,
            oDir = "mapping_qc/reports"
        }
    }
    if(runVerifyBamIdReport){
        call verifyBamIdReport {
            input:
                verifyBamIdSelfSM = verifyBamIdSelfSM,
                docker = gencompassDocker,
                project = project,
                # date = extractDate.date,
                oDir = "mapping_qc/reports"
        }
    }
    
    if(runKrakenReport){
        call krakenReport {
        input:
            krakenReports = read_lines(select_first([krakenReportFOF, buildFileList.krakenReportFiles])),
            krakenOutputs = read_lines(select_first([krakenOutputFOF, buildFileList.krakenOutputFiles])),
            docker = gencompassDocker,
            project = project,
            # date = extractDate.date, 
            oDir = "mapping_qc/reports"
    }

    }
    

    output {
        File? somalierRelateSamplesTSV = somalierRelate.samplesTSV
        File? somalierRelateSampleHTML = somalierRelate.sampleHTML
        File? somalierRelateRelatednessPairs = somalierRelate.relatednessPairs
        File? baitBiasReport = sequencingArtifactReport.baitBiasReport
        File? preAdapterReport = sequencingArtifactReport.preAdapterReport
        File? mappingParabricksReportXLSX = mappingParabricksReport.mappingReport
        File? multiqcHTML = multiqcReport.multiqcHTML
        File? multiqcDataJson= multiqcReport.multiqcDataJson
        File? verifyBamIdReportXLSX = verifyBamIdReport.verifyBamIdReportXLSX 
        File? somalierRelatednessReport = somalierSexAndRelatednessReport.somalierRelatednessReport
        File? somalierSexReport = somalierSexAndRelatednessReport.somalierSexReport
        File? krakenReportXLSX = krakenReport.krakenReportXLSX
        
    }

}

task extractDate {
    command <<<
        date=$(date '+%Y%m%d')
        echo $date
    >>>
    output {
        String date = read_lines(stdout())[0]
    }
}

task buildFileList {
    input{
        File sampleList
        String workflowResultsDirectory
        String docker
    }
    command <<<
    while read sample; do 
        echo "~{workflowResultsDirectory}/mapping_qc/bammetrics/$sample.bammetrics.txt" >> bammetrics_files.txt
        echo "~{workflowResultsDirectory}/mapping_qc/collectmultiplemetrics/$sample/sequencingArtifact.bait_bias_summary_metrics.txt" >> cmm_bait_bias_files.txt
        echo "~{workflowResultsDirectory}/mapping_qc/collectmultiplemetrics/$sample/sequencingArtifact.pre_adapter_summary_metrics.txt" >> cmm_pre_adapter_files.txt

        echo "~{workflowResultsDirectory}/mapping_qc/collectmultiplemetrics/$sample/alignment.txt" >> cmm_summary_files.txt
        echo "~{workflowResultsDirectory}/mapping_qc/collectmultiplemetrics/$sample/gcbias_detail.txt" >> cmm_summary_files.txt
        echo "~{workflowResultsDirectory}/mapping_qc/collectmultiplemetrics/$sample/gcbias_summary.txt" >> cmm_summary_files.txt
        echo "~{workflowResultsDirectory}/mapping_qc/collectmultiplemetrics/$sample/insert_size.txt" >> cmm_summary_files.txt
        echo "~{workflowResultsDirectory}/mapping_qc/collectmultiplemetrics/$sample/mean_quality_by_cycle.txt" >> cmm_summary_files.txt
        echo "~{workflowResultsDirectory}/mapping_qc/collectmultiplemetrics/$sample/qualityscore.txt" >> cmm_summary_files.txt
        echo "~{workflowResultsDirectory}/mapping_qc/collectmultiplemetrics/$sample/quality_yield.txt" >> cmm_summary_files.txt

        echo "~{workflowResultsDirectory}/mapping_qc/kraken2/$sample/$sample.kraken.output" >> kraken_output_files.txt
        echo "~{workflowResultsDirectory}/mapping_qc/kraken2/$sample/$sample.kraken.mpaStyle.report" >> kraken_report_files.txt

        echo "~{workflowResultsDirectory}/mapping_qc/somalier/extract/$sample.somalier" >> somalier_extract_files.txt

        echo "~{workflowResultsDirectory}/mapping_qc/verifybamid/$sample/$sample.selfSM" >> verifybamid_selfsm_files.txt
    done < ~{sampleList}
    >>>

    output {
        File bammetricsFiles = "bammetrics_files.txt"
        File cmmBaitBiasFiles = "cmm_bait_bias_files.txt"
        File cmmPreAdapterFiles = "cmm_pre_adapter_files.txt"
        File cmmSummaryFiles = "cmm_summary_files.txt"
        File krakenOutputFiles = "kraken_output_files.txt"
        File krakenReportFiles = "kraken_report_files.txt"
        File somalierExtractFiles = "somalier_extract_files.txt"
        File verifybamidSelfSmFiles = "verifybamid_selfsm_files.txt"


    }
    runtime {
        docker : docker
        disks : "local-disk 1 SSD"
        cpu : 2
        memory : "1 GB"
        hpcMemory : 1
        hpcQueue : "norm"
        hpcRuntimeMinutes : 30
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : 3
    }



}

task somalierRelate {
    input {
        Array[File] extracted
        String project
        
        String docker ="cgrlab/somalier@sha256:3f63c21f6741e49a13dca4384840d68802eec8ead950771629bc233a055a6a1e"
        Int nThreads=8
        String hpcQueue="norm"
        Int memoryGiB=2
        Int diskGB=0
        Int runtimeMinutes=60
        Int maxPreemptAttempts=3

    }
    String oDir = "mapping_qc/somalier/relate"
    Int auto_diskGB = if diskGB < 1 then 500 else diskGB
    # String localTarball = basename(labeledSamplesTarball)
    # String labeledSamplesDir = basename(labeledSamplesTarball, ".tar")
    # mv ~{labeledSamplesTarball} ${localTarball}
    # tar xvf ~{localTarball}
    command <<<
        set -euo pipefail
        mkdir -p ~{oDir} 
        cd ~{oDir}

        somalier relate --output-prefix ~{project} ~{sep=" " extracted} 
    >>>

    output {
        File samplesTSV = "~{oDir}/~{project}.samples.tsv"
        File sampleHTML = "~{oDir}/~{project}.html"
        File relatednessPairs = "~{oDir}/~{project}.pairs.tsv"
        
    }
    runtime {
        docker : docker
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{memoryGiB} GB"
        hpcMemory : memoryGiB
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

#-------------------------------------
# Sequencing Artifact Reports
#   - pre
#   - bait
#-------------------------------------
task sequencingArtifactReport {
    input {
        Array[File] preAdapterSummaryMetrics 
        Array[File] baitBiasSummaryMetrics
        String project
        String oDir
        String docker
        # String date
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        String hpcQueue = "norm"
        Int memoryGiB = 20

    }
    Int autoDiskGB = ceil(size(preAdapterSummaryMetrics, "GiB")) + ceil(size(baitBiasSummaryMetrics, "GiB"))
    command <<<
    set -euo pipefail
    mkdir -p ~{oDir} 
    sequencing_artifact_reports.py \
        --seqart_pre_list ~{write_lines(preAdapterSummaryMetrics)} \
        --seqart_bait_list ~{write_lines(baitBiasSummaryMetrics)} \
        --project ~{project} \
        -odir ~{oDir}
    
    mv ~{oDir}/~{project}.bam_QC_C_sequencingArtifact.bait_bias-*.xlsx ~{oDir}/~{project}.bam_QC_C_sequencingArtifact.bait_bias.xlsx
    mv ~{oDir}/~{project}.bam_QC_C_sequencingArtifact.pre_adapter-*.xlsx ~{oDir}/~{project}.bam_QC_C_sequencingArtifact.pre_adapter.xlsx
    >>>
    
    output {
        File baitBiasReport = "~{oDir}/~{project}.bam_QC_C_sequencingArtifact.bait_bias.xlsx"
        File preAdapterReport = "~{oDir}/~{project}.bam_QC_C_sequencingArtifact.pre_adapter.xlsx"
        # Array[File] sequencingArtifactReportsXLSX = glob("~{oDir}/*.xlsx")
    }
    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} SSD"
        cpu : 2
        memory : "~{memoryGiB} GiB"
        hpcMemory : memoryGiB
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

#-------------------------------------
# BamMetrics Report
#-------------------------------------
task bammetricsReport {
    input {
        Array[File] bammetricsIndividualReports 
        String project
        String oDir
        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        String hpcQueue = "norm"
        Int memoryGiB = 20

    }
    Int autoDiskGB = ceil(size(bammetricsIndividualReports, "GiB") * 1.2)
    command <<<
        set -euo pipefail
        mkdir -p ~{oDir} 
        bammetrics_report.py \
        --bammetrics_list ~{write_lines(bammetricsIndividualReports)} \
        --project ~{project} \
        -odir ~{oDir}
    >>>
    
    output {
        File bammetricsReportTable = "~{oDir}/~{project}.bammetrics_report.tsv"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} SSD"
        cpu : 2
        memory : "~{memoryGiB} GiB"
        hpcMemory : memoryGiB
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}
#-------------------------------------
# MultiQC report
#-------------------------------------
task multiqcReport {
    input {
        Array[File] inputData 
        String filename
        String title
        String? comment 
        String oDir
        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        String hpcQueue = "norm"
        Int memoryGiB = 20

    }
    Int autoDiskGB = ceil(size(inputData, "GiB") * 1.2)
    command <<<
            set -euo pipefail
            mkdir -p ~{oDir} 
            multiqc \
            --outdir ~{oDir} \
            --filename ~{filename} \
            --title ~{title} \
            ~{"--comment " + comment} \
            --file-list ~{write_lines(inputData)}
    >>>
    
    output {
        File multiqcHTML = "~{oDir}/~{filename}.html"
        File multiqcDataJson="~{oDir}/~{filename}_data/multiqc_data.json"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} SSD"
        cpu : 2
        memory : "~{memoryGiB} GiB"
        hpcMemory : memoryGiB
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}


#-------------------------------------
# Post Mapping Parabricks report
#-------------------------------------
task mappingParabricksReport {
    input {
        File multiqcJSON
        File bammetricsReportTable
        String project
        String oDir
        String docker
        # String date
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        String hpcQueue = "norm"
        Int memoryGiB = 20
    }
    Int autoDiskGB = ceil(size(multiqcJSON, "GiB")*2) + ceil(size(bammetricsReportTable, "GiB")*2)
    command <<<
        set -euo pipefail
        mkdir -p ~{oDir} 
        post_mapping_report.py \
            --cmm_multiqc_json ~{multiqcJSON} \
            --bammetrics_report ~{bammetricsReportTable} \
            --project ~{project} \
            -odir ~{oDir}
        
        mv ~{oDir}/~{project}.post_mapping_report-*.xlsx ~{oDir}/~{project}.post_mapping_report.xlsx

    >>>
    
    output {
        File mappingReport = "~{oDir}/~{project}.post_mapping_report.xlsx"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} SSD"
        cpu : 2
        memory : "~{memoryGiB} GiB"
        hpcMemory : memoryGiB
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}
#-------------------------------------
# VerifyBamID2 Report
#-------------------------------------
task verifyBamIdReport {
    input {
        Array[File] verifyBamIdSelfSM
        String project
        String oDir
        String docker
        # String date
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        String hpcQueue = "norm"
        Int memoryGiB = 20
    }
    Int autoDiskGB = ceil(size(verifyBamIdSelfSM, "GiB") * 1.2) 
    command <<<
        set -euo pipefail
        mkdir -p ~{oDir} 
        verifybamid2_report.py \
            --selfsm_list ~{write_lines(verifyBamIdSelfSM)} \
            --project ~{project} \
            -odir ~{oDir}
        
        mv ~{oDir}/~{project}.verifybamid2.selfSM-*.xlsx ~{oDir}/~{project}.verifybamid2.selfSM.xlsx

    >>>
    
    output {
        File verifyBamIdReportXLSX = "~{oDir}/~{project}.verifybamid2.selfSM.xlsx"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} SSD"
        cpu : 2
        memory : "~{memoryGiB} GiB"
        hpcMemory : memoryGiB
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

#-------------------------------------
# Somalier Sex and Relatedness Reports
#-------------------------------------
task somalierSexAndRelatednessReport {
    input {
        File somalierSamples
        File somalierPairs
        Int ibs2Cutoff = 5000
        Float relatednessCutoff = 0.0
        String project
        String oDir
        # String date
        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        String hpcQueue = "norm"
        Int memoryGiB = 20
    }
    
    Int autoDiskGB = ceil(size(somalierSamples, "GiB")*1.2 + size(somalierPairs, "GiB")*1.2)
    command <<<
        set -euo pipefail
        mkdir -p ~{oDir} 
        somalier_report.py \
            --somalier_samples ~{somalierSamples} \
            --somalier_pairs ~{somalierPairs} \
            --ibs2_cutoff ~{ibs2Cutoff} \
            --relatedness_cutoff ~{relatednessCutoff} \
            --project ~{project} \
            -odir ~{oDir}
        
        mv ~{oDir}/~{project}.bam_QC_somalier_relatedness-*.xlsx ~{oDir}/~{project}.bam_QC_somalier_relatedness.xlsx
        mv ~{oDir}/~{project}.bam_QC_somalier_sex_report-*.xlsx ~{oDir}/~{project}.bam_QC_somalier_sex_report.xlsx
    >>>
    
    output {
        File somalierRelatednessReport = "~{oDir}/~{project}.bam_QC_somalier_relatedness.xlsx"
        File somalierSexReport = "~{oDir}/~{project}.bam_QC_somalier_sex_report.xlsx"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} SSD"
        cpu : 2
        memory : "~{memoryGiB} GiB"
        hpcMemory : memoryGiB
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}


#-------------------------------------
# Kraken2 Report
#-------------------------------------
task krakenReport {
    input {
        Array[File] krakenReports
        Array[File] krakenOutputs
        String project

        String oDir
        String docker
        # String date
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        String hpcQueue = "norm"
        Int memoryGiB = 20
    }
    Int autoDiskGB = ceil(size(krakenReports, "GiB")*1.2 + size(krakenOutputs, "GiB")*1.2)
    command <<<
        set -euo pipefail
        mkdir -p ~{oDir} 
        kraken_report.py \
            --kraken_reports ~{write_lines(krakenReports)} \
            --kraken_outputs ~{write_lines(krakenOutputs)} \
            --project ~{project} \
            -odir ~{oDir}
        
        mv ~{oDir}/~{project}.bam-QC-Kraken-genus-level-*.xlsx ~{oDir}/~{project}.bam-QC-Kraken-genus-level.xlsx
    >>>
    
    output {
        File krakenReportXLSX = "~{oDir}/~{project}.bam-QC-Kraken-genus-level.xlsx"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} SSD"
        cpu : 2
        memory : "~{memoryGiB} GiB"
        hpcMemory : memoryGiB
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}
