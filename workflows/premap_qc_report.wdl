version 1.0

workflow PremapReport {
    input {

        Boolean runFastpReport = true
        Boolean runFastqScreenReport = true
        Boolean runMultiQCReport = true

        Boolean includeFastpInMultiQC = true
        Boolean includeFastqScreenInMultiQC = true
        Boolean includeFastqcInMultiQC = true

        File? sampleFastqLocations
        File? sampleList
        String? workflowResultsDirectory

        File? fastpJsonFOF
        File? fastqScreenTxtFOF
        File? fastQcZipFOF

        File manifest
        String project
        String sequencingType="Illumina"

        String gencompassDocker
        String multiqcDocker
        RuntimeAttributes? fastpReportRuntimeAttributes
        RuntimeAttributes? fastqScreenReportRuntimeAttributes
        RuntimeAttributes? multiqcReportRuntimeAttributes
    }
    if(defined(sampleFastqLocations) && defined(sampleList) && defined(workflowResultsDirectory)){
        call buildFileLists {
            input:
                sampleFastqLocations=select_first([sampleFastqLocations]),
                manifest=manifest,
                sampleList=select_first([sampleList]),
                workflowResultsDirectory=select_first([workflowResultsDirectory]),
                docker=gencompassDocker
        }
    }
    if( runFastpReport ||(runMultiQCReport && includeFastpInMultiQC)){
        Array[File] fastpJsonArray = read_lines(select_first([fastpJsonFOF, buildFileLists.fastpLocations]))
    }
    if( runFastqScreenReport ||(runMultiQCReport && includeFastqScreenInMultiQC)){
        Array[File] fastqScreenTxtArray = read_lines(select_first([fastqScreenTxtFOF, buildFileLists.fastqScreenLocations]))
    }
    if(runMultiQCReport && includeFastqcInMultiQC){
        Array[File] fastqcZipArray = read_lines(select_first([fastQcZipFOF, buildFileLists.fastqcLocations]))
    }
    

    if(runFastpReport){
        call fastpReport {
            input:
                fastpJsonArray=select_first([fastpJsonArray]),
                manifest=manifest,
                sequencingType=sequencingType,
                projectName=project,
                oDir="premap_qc/reports",
                docker=gencompassDocker,
                runtimeAttributes=fastpReportRuntimeAttributes
        }

    }
    
    if(runFastqScreenReport){
        call fastqScreenReport {
            input:
                fastqScreenTxtArray=select_first([fastqScreenTxtArray]),
                projectName=project,
                oDir="premap_qc/reports",
                docker=gencompassDocker,
                runtimeAttributes=fastqScreenReportRuntimeAttributes
        }
    }
    
    if(runMultiQCReport){
        call multiqcReport {
            input:
                title="~{project}-PremapQC-MultiQC-Report",
                docker = multiqcDocker,
                inputData = flatten([select_first([fastpJsonArray]),select_first([fastqScreenTxtArray]), select_first([fastqcZipArray])]),
                oDir="premap_qc/reports",
                filename = "premap_qc_multiqc",
                runtimeAttributes=multiqcReportRuntimeAttributes
        }
    }
    output{
        File? fastpReportFile = fastpReport.fastpReportFile
        File? fastqScreenReportFile = fastqScreenReport.fastqScreenReportFile
        File? multiqcHTML = multiqcReport.multiqcHTML
        File? multiqcDataJson = multiqcReport.multiqcDataJson
    }

}

task buildFileLists {
    input {
        File sampleList
        File sampleFastqLocations
        File manifest
        String workflowResultsDirectory
        String docker
    }
    command<<<
    set -euxo pipefail
    create_premap_report_file_list.py \
    --sample_list ~{sampleList} \
    --manifest ~{manifest} \
    --sample_fastq_files ~{sampleFastqLocations} \
    --workflow_results_dir ~{workflowResultsDirectory}
    >>>

    output {
        File fastpLocations = "fastp_file_list.txt"
        File fastqcLocations = "fastqc_files.txt"
        File fastqScreenLocations = "fastq_screen_files.txt"
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

task extractDate {
    command <<<
        date=$(date '+%Y%m%d')
        echo $date
    >>>
    output {
        String date = read_lines(stdout())[0]
    }
}

task fastpReport {
    input{
        Array[File] fastpJsonArray
        File manifest
        String projectName
        String sequencingType
        String oDir
        # String date
        String docker
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 10,
                    "cpuCount" : 2,
                    "diskGiB" : 0,
                    "runtimeMinutes": 240,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(size(fastpJsonArray, "GiB") * 1.5 )  else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])

    command <<<
        set -euxo pipefail
        mkdir -p ~{oDir}
        fastp_report.py \
        --fastp_json_list ~{write_lines(fastpJsonArray)} \
        --manifest ~{manifest} \
        --project_name ~{projectName} \
        --sequencing_type ~{sequencingType} \
        --output_directory ~{oDir}

        mv ~{oDir}/~{projectName}.pre-mapping-QC-A-fastp_report-*.xlsx ~{oDir}/~{projectName}.pre-mapping-QC-A-fastp_report.xlsx

    >>>

    output {
        File fastpReportFile = "~{oDir}/~{projectName}.pre-mapping-QC-A-fastp_report.xlsx"
    }

    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} ~{select_first([runtimeAttributesOverride.diskType, defaultRuntimeAttributes.diskType])}"
        cpu : select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
        memory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) + " GiB"
        
        hpcMemory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB])
        memory_mb : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) * 1024
        hpcQueue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        queue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        hpcRuntimeMinutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])
        runtime_minutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : select_first([runtimeAttributesOverride.maxPreemptAttempts, defaultRuntimeAttributes.maxPreemptAttempts])
    }
}

task fastqScreenReport {
    input{
        Array[File] fastqScreenTxtArray
        String projectName
        String oDir
        # String date

        String docker
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 10,
                    "cpuCount" : 2,
                    "diskGiB" : 0,
                    "runtimeMinutes": 120,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(size(fastqScreenTxtArray, "GiB") * 1.5 ) else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])

    command <<<
        set -euxo pipefail
        mkdir -p ~{oDir} 
        fastq_screen_report.py \
        --fastq_screen_list ~{write_lines(fastqScreenTxtArray)} \
        --project_name ~{projectName} \
        --output_directory ~{oDir}

        mv ~{oDir}/~{projectName}.pre-mapping-QC-B-fastq_screen_report*.xlsx ~{oDir}/~{projectName}.pre-mapping-QC-B-fastq_screen_report.xlsx 
    >>>

    output {
        File fastqScreenReportFile = "~{oDir}/~{projectName}.pre-mapping-QC-B-fastq_screen_report.xlsx"
        
    }

    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} ~{select_first([runtimeAttributesOverride.diskType, defaultRuntimeAttributes.diskType])}"
        cpu : select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
        memory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) + " GiB"
        
        hpcMemory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB])
        memory_mb : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) * 1024
        hpcQueue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        queue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        hpcRuntimeMinutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])
        runtime_minutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : select_first([runtimeAttributesOverride.maxPreemptAttempts, defaultRuntimeAttributes.maxPreemptAttempts])
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
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 100,
                    "cpuCount" : 2,
                    "diskGiB" : 0,
                    "runtimeMinutes": 240,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int autoDiskGB = ceil(size(inputData, "GiB") * 2)
    command <<<
        set -euxo pipefail
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
        disks : "local-disk ~{autoDiskGB} ~{select_first([runtimeAttributesOverride.diskType, defaultRuntimeAttributes.diskType])}"
        cpu : select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
        memory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) + " GiB"
        
        hpcMemory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB])
        memory_mb : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) * 1024
        hpcQueue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        queue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        hpcRuntimeMinutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])
        runtime_minutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : select_first([runtimeAttributesOverride.maxPreemptAttempts, defaultRuntimeAttributes.maxPreemptAttempts])
    }
}

struct RuntimeAttributes {
    Int? memoryGiB
    Int? cpuCount
    Int? diskGiB
    Int? bootDiskGiB
    Int? preemptibleTries
    Int? maxRetries
    String? acceleratorType
    Int? acceleratorCount
    String? acceleratorDriverVersion
    Int? maxPreemptAttempts
    Int? runtimeMinutes
    String? hpcQueue
    String? diskType
}