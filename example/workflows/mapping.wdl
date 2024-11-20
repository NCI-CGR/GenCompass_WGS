version 1.0

workflow Mapping {
    # Maps fastq files to bam files and performs QC analysis. 
    # If a mapped bam file is provided instead QC will be 
    # performed on provided bam file. 
    input {
       
        String sampleID
        File referenceTarball 

        # Array[Array[File]] pairedReads
        File manifest
        File? sampleFastqLocations
        Array[File] sampleFastqFiles=[]
        File? mappedBAM
        File? mappedBAI

        Boolean runBamMetrics=true
        Boolean runCollectMultipleMetrics=true
        Boolean runKraken=true
        Boolean runSamtoolsCoverage=true
        Boolean runSomalier=true
        Boolean runVerifyBamID=true

        File? fq2bamKnownSites
        Boolean fq2bamUseBestPractices=false
        String fq2bamToolOptions="--gpusort --gpuwrite --low-memory"

        Float kraken2Confidence=0.2
        File? kraken2DatabaseTarball

        File? somalierExtractSites

        File? verifybamidSvdTarball

        
        String gencompassDocker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        String parabricksDocker="nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"
        String? verifybamidDocker
        String? krakenDocker
        String? somalierDocker

        RuntimeAttributes fq2bamRuntimeAttributes = {
                    "memoryGiB" : 192,
                    "cpuCount" : 48,
                    "acceleratorType" : "nvidia-tesla-a10g",
                    "acceleratorCount": 4,
                    "diskGiB" : 0
        }

        RuntimeAttributes collectmultiplemetricsRuntimeAttributes = {
                    "memoryGiB" : 32,
                    "cpuCount" : 8,
                    "acceleratorType" : "nvidia-tesla-a10g",
                    "acceleratorCount": 1,
                    "diskGiB" : 0
        }

        RuntimeAttributes bammetricsRuntimeAttributes = {
                    "memoryGiB" : 64,
                    "cpuCount" : 16,
                    "acceleratorType" : "nvidia-tesla-a10g",
                    "acceleratorCount": 1,
                    "diskGiB" : 0
        }

        RuntimeAttributes verifybamidRuntimeAttributes = {
                    "memoryGiB" : 16,
                    "cpuCount" : 8,
                    "diskGiB" : 0
        }
        RuntimeAttributes samtoolsCoverageRuntimeAttributes = {
                    "memoryGiB" : 8,
                    "cpuCount" : 2,
                    "diskGiB" : 0
        }

        RuntimeAttributes extractUnmappedReadsRuntimeAttributes = {
                    "memoryGiB" : 8,
                    "cpuCount" : 2,
                    "diskGiB" : 0
        }

        RuntimeAttributes krakenRuntimeAttributes = {
                    "memoryGiB" : 16,
                    "cpuCount" : 4,
                    "diskGiB" : 0
        }
        RuntimeAttributes somalierExtractRuntimeAttributes = {
                    "memoryGiB" : 4,
                    "cpuCount" : 2,
                    "diskGiB" : 0
        }
    }
    if( defined(sampleFastqLocations) || length(sampleFastqFiles) > 0 )
    {
        
        # File fastqFiles = select_first([read_lines(sampleFastqLocations), write_lines(sampleFastqFiles)])
        # Array[String] fastqFiles = read_lines(select_first([sampleFastqLocations, write_lines(sampleFastqFiles)])) 

        File fastqFileList = select_first([sampleFastqLocations, write_lines(sampleFastqFiles)])
        # Array[File] fastqFilesFromFOF = read_lines(select_first([sampleFastqLocations, []]))
        # Array[File] fastqFiles = flatten([fastqFilesFromFOF, sampleFastqFiles])


        call prepare_fq2bam_list{
            input:
                fastqFileList=fastqFileList,
                manifest=manifest,
                sampleID=sampleID,
                docker=gencompassDocker
        }
        Array[File] preparedSampleFastqFiles = read_lines(prepare_fq2bam_list.sampleFastqFiles)

        if(defined(fq2bamKnownSites)){
            File fq2bamKnownSitesTBI=select_first([fq2bamKnownSites])+".tbi"
        }
        
        call multifq2bam as fq2bam {
            input:
                fastqReadGroupList=prepare_fq2bam_list.fq2bamList,
                pairedReads=preparedSampleFastqFiles,
                sampleID=sampleID,
                inputRefTarball=referenceTarball,
                inputKnownSitesVCF=fq2bamKnownSites,
                inputKnownSitesTBI=fq2bamKnownSitesTBI,
                useBestPractices=fq2bamUseBestPractices,
                fq2bamToolOptions=fq2bamToolOptions,
                docker=parabricksDocker,
                runtimeAttributes = fq2bamRuntimeAttributes
        }
    }
    if(defined(mappedBAM)){
        File autoBAI = select_first([mappedBAI, "~{mappedBAM}.bai"])
        BamFile mappedBamProvided={
            "bam" : select_first([mappedBAM]),
            "bai": autoBAI
        }
        
    }
    BamFile mappedBamQC = select_first([fq2bam.bamFile, mappedBamProvided ])
    # Array[File] fastqFiles = read_lines(sampleFastqLocations)

    if(runCollectMultipleMetrics){
        call collectmultiplemetrics{
            input:
                inputBAM=mappedBamQC.bam,
                inputRefTarball=referenceTarball,
                sampleID=sampleID,
                docker=parabricksDocker,
                runtimeAttributes = collectmultiplemetricsRuntimeAttributes
        }

    }
    if(runBamMetrics){

        call bammetrics{
            input:
                inputBAM=mappedBamQC.bam,
                inputRefTarball=referenceTarball,
                sampleID=sampleID,
                docker=parabricksDocker,
                runtimeAttributes=bammetricsRuntimeAttributes
        }
    }
    
    if(runVerifyBamID){
        call verifybamid {
            input:
                sampleID=sampleID,
                bamFile=mappedBamQC,
                referenceTarball=referenceTarball,
                svdTarball=select_first([verifybamidSvdTarball]),
                docker=verifybamidDocker,
                runtimeAttributes=verifybamidRuntimeAttributes
        }
    }
    
    if(runSamtoolsCoverage){
        call samtools_coverage {
            input:
                sampleID=sampleID,
                bamFile=mappedBamQC,
                docker=gencompassDocker,
                runtimeAttributes=samtoolsCoverageRuntimeAttributes
                # bamIndex=fq2bam.outputBAI, 
        } 
    }
    
    if(runKraken){
        call extract_unmapped_reads {
            input:
                sampleID=sampleID,
                bam=mappedBamQC.bam,
                docker=gencompassDocker,
                runtimeAttributes=extractUnmappedReadsRuntimeAttributes
        }
        call kraken2{
            
            input:
                sampleID=sampleID,
                unmappedReads=extract_unmapped_reads.unmappedReads,
                oDir="mapping_qc/kraken2/~{sampleID}",
                databaseTarball=select_first([kraken2DatabaseTarball]),
                confidence=kraken2Confidence,
                docker=krakenDocker,
                runtimeAttributes=krakenRuntimeAttributes
        }
    }
    if(runSomalier){
        call somalier_extract {
            input:
                bamFile=mappedBamQC,
                referenceTarball=referenceTarball,
                sites=select_first([somalierExtractSites]),
                docker=somalierDocker,
                runtimeAttributes=somalierExtractRuntimeAttributes
        }

    }
    
    output{
        BamFile? fq2bamBAM=fq2bam.bamFile
        # File fq2bamBAI=fq2bam.outputBAI
        File? fq2bamBQSR=fq2bam.outputBQSR
        File? fq2bamMetrics = fq2bam.metrics
        Array[File]? multipleMetrics = collectmultiplemetrics.multipleMetrics
        File? bammetricsFile = bammetrics.oMetrics
        File? coverage = samtools_coverage.coverage

        File? selfSM = verifybamid.selfSM
        File? ancestory = verifybamid.ancestory

        File? krakenReport = kraken2.krakenReport
        File? krakenOutput = kraken2.krakenOutput
        File? somalierExtract = somalier_extract.extracted
    }
}

# PREPARE FQ2BAM LIST WITH READ GROUP INFO
task prepare_fq2bam_list {
    input{
        File fastqFileList
        File manifest
        String sampleID
        String docker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        
        String hpcQueue="norm"
        Int gbRAM=2
        Int diskGB=0
        Int runtimeMinutes=15
        Int maxPreemptAttempts=3
    }

    String oDir = "fq2bam/~{sampleID}"
    
    Int autoDiskGB = if diskGB < 1 then ceil(2.0 * size(manifest,  "GiB"))   + 5 else diskGB

    command {
        set -euxo pipefail
        mkdir -p ~{oDir}
        
        prepare_fq2bam_list.py \
        --manifest ~{manifest} \
        --fastq_files ~{fastqFileList} \
        --sample ~{sampleID} \
        --output_directory ~{oDir}
    }
    output {
        File fq2bamList = "~{oDir}/~{sampleID}.fq2bam_list.txt"
        File sampleFastqFiles = "~{oDir}/~{sampleID}_fastq_files.txt"
    }
    runtime {
        docker: docker
        disks : "local-disk ~{autoDiskGB} SSD"
        cpu : 2
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

# PARABRICKS FQ2BAM
task multifq2bam {
    input {

        Array[File] pairedReads
        File fastqReadGroupList

        File inputRefTarball
        String sampleID = "SAMPLE"

        File? inputKnownSitesVCF
        File? inputKnownSitesTBI
        String docker = "nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"
        Boolean useBestPractices = false
        String fq2bamToolOptions="--gpusort --gpuwrite --low-memory"

        RuntimeAttributes? runtimeAttributes

        String tmpDir = "tmp_fq2bam"
    }

    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 192,
                    "cpuCount" : 48,
                    "acceleratorType" : "nvidia-tesla-a10g",
                    "acceleratorCount": 4,
                    "diskGiB" : 0,
                    "diskType" : "SSD",
                    "runtimeMinutes": 600,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "gpu",
                    "acceleratorDriverVersion": "525.60.13"
        }
    
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(4.0 * size(pairedReads,  "GB")) + ceil(4.0 * size(inputRefTarball,  "GB")) + ceil(4.0 * size(inputKnownSitesVCF,  "GB")) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 
    String best_practice_args = if useBestPractices then "--bwa-options \" -Y -K 100000000 \" " else ""
    
    String ref = basename(inputRefTarball, ".tar")
    String oDir = "fq2bam/~{sampleID}"

    String bqsrFilename = "~{oDir}/~{sampleID}.BQSR-REPORT.txt"
    String bamFilename = "~{oDir}/~{sampleID}.bam"
    String metricsFilename="~{oDir}/~{sampleID}.metrics.txt"
    


    command {
        set -euxo pipefail

        mkdir -p ~{oDir}
        cp ~{sep=" " pairedReads} ./

        mkdir -p ~{tmpDir} && \
        tar xf ~{inputRefTarball} && \
        pbrun fq2bam \
            --tmp-dir ~{tmpDir} \
            --in-fq-list ~{fastqReadGroupList} \
            --ref ~{ref} \
            ~{"--knownSites " + inputKnownSitesVCF + " --out-recal-file " + bqsrFilename} \
            ~{best_practice_args} \
            --out-bam ~{bamFilename} \
            --out-duplicate-metrics ~{metricsFilename}  \
            ~{fq2bamToolOptions}

        rm ~{ref}*
        rm *.fastq.gz
    }

    output {
        BamFile bamFile = {
            "bam" : "~{bamFilename}",
            "bai": "~{bamFilename}.bai"
        }
        # File outputBAM = "~{bamFilename}"
        # File outputBAI = "~{bamFilename}.bai"
        File? outputBQSR = "~{bqsrFilename}"
        File metrics = "~{metricsFilename}"
        
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

        gpuType : select_first([runtimeAttributesOverride.acceleratorType, defaultRuntimeAttributes.acceleratorType])
        gpuCount : select_first([runtimeAttributesOverride.acceleratorCount, defaultRuntimeAttributes.acceleratorCount])

        acceleratorCount: select_first([runtimeAttributesOverride.acceleratorCount, defaultRuntimeAttributes.acceleratorCount])
        acceleratorType: select_first([runtimeAttributesOverride.acceleratorType, defaultRuntimeAttributes.acceleratorType])

        nvidiaDriverVersion : select_first([runtimeAttributesOverride.acceleratorDriverVersion, defaultRuntimeAttributes.acceleratorDriverVersion])
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : select_first([runtimeAttributesOverride.maxPreemptAttempts, defaultRuntimeAttributes.maxPreemptAttempts])
    }
}

# KRAKEN EXTRACT UNMAPPED READS
task extract_unmapped_reads {
    input {
        File bam
        String sampleID
        
        String docker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 50,
                    "cpuCount" : 2,
                    "diskGiB" : 0,
                    "runtimeMinutes": 60,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(2.0 * size(bam,  "GB")) + 25 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    command {
        set -euxo pipefail
        # Extract unmapped reads (samtools -f 4) - note this breaks up the read pairing
        # Run samtools fastq to generate fastq from unmapped bam
        samtools view -u -f 4 ~{bam} | samtools fastq >  ~{sampleID}.unmapped.fq

    }
    output {
        File unmappedReads = "~{sampleID}.unmapped.fq"
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
# KRAKEN2
task kraken2 {
    input {
        File databaseTarball
        File unmappedReads
        String sampleID

        Float confidence = 0.2
        String oDir = "mapping_qc/kraken2"

        String docker="cgrlab/kraken@sha256:63d241d3819fe82157abbac68472e6a8bbb6b02d7f5175845175887de554bb98"

        RuntimeAttributes? runtimeAttributes
    }
    
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 50,
                    "cpuCount" : 8,
                    "diskGiB" : 0,
                    "runtimeMinutes": 60,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])< 1 then ceil(2.0 * size(databaseTarball,  "GB")) + ceil(2.0 * size(unmappedReads,  "GB")) + 25 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])

    String localTarball = basename(databaseTarball)
    
    String dbName = basename(basename(databaseTarball, ".tar.gz"), ".tar")

    command {
        set -euxo pipefail
        mkdir -p ~{oDir}
        ln -s ~{databaseTarball} ${localTarball} && \
        tar xvf ~{localTarball}

        
        # Run Kraken
        kraken2 --db ~{dbName} \
                    --confidence ~{confidence} \
                    --threads ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} \
                    --output ~{oDir}/~{sampleID}.kraken.output \
                    --report ~{oDir}/~{sampleID}.kraken.mpaStyle.report \
                    --use-mpa-style \
                    ~{unmappedReads}
        
        rm -r ~{dbName}
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
    output{
        File krakenOutput="~{oDir}/~{sampleID}.kraken.output"
        File krakenReport="~{oDir}/~{sampleID}.kraken.mpaStyle.report"
    }
}

# PARABRICKS COLLECTMULTIPLEMETRICS
task collectmultiplemetrics {
    input {
        File inputBAM
        File inputRefTarball

        String sampleID = ""

        String docker = "nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"


        RuntimeAttributes? runtimeAttributes 
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 64,
                    "cpuCount" : 16,
                    "acceleratorType" : "nvidia-tesla-a10g",
                    "acceleratorCount": 1,
                    "diskGiB" : 0,
                    "runtimeMinutes": 240,
                    "maxPreemptAttempts": 3,
                    "diskType" : "SSD",
                    "hpcQueue": "gpu",
                    "acceleratorDriverVersion": "525.60.13"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])

    # String outbase = basename(inputBAM, ".bam")
    String autoSampleID = if sampleID == "" then basename(basename(inputBAM, ".bam"), ".cram") else sampleID
    String localTarball = basename(inputRefTarball)
    String ref = basename(inputRefTarball, ".tar")

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(2.0 * size(inputBAM,  "GB")) +ceil(2.0 * size(inputRefTarball,  "GB")) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 

    String oDir = "mapping_qc/collectmultiplemetrics/~{autoSampleID}" 

    command<<<
        # mkdir -p ~{oDir}
        set -euxo pipefail
        ln -s ~{inputRefTarball} ~{localTarball} && \
        tar xvf ~{localTarball}

        pbrun collectmultiplemetrics \
        --ref ~{ref} \
        --bam ~{inputBAM} \
        --out-qc-metrics-dir ~{oDir} \
        --gen-all-metrics

        rm ~{ref}*

        ls ~{oDir} | awk -v prefix="~{oDir}/" '$0=prefix$0' > files.txt
    
    >>>

    output {

       
        Array[File] multipleMetrics = [
             "~{oDir}/alignment.txt",
            "~{oDir}/base_distribution_by_cycle.pdf",
            "~{oDir}/base_distribution_by_cycle.png",
            "~{oDir}/base_distribution_by_cycle.txt",
            "~{oDir}/gcbias.pdf",
            "~{oDir}/gcbias_0.png",
            "~{oDir}/gcbias_detail.txt",
            "~{oDir}/gcbias_summary.txt",
            "~{oDir}/insert_size.pdf",
            "~{oDir}/insert_size.png",
            "~{oDir}/insert_size.txt",
            "~{oDir}/mean_quality_by_cycle.pdf",
            "~{oDir}/mean_quality_by_cycle.png",
            "~{oDir}/mean_quality_by_cycle.txt",
            "~{oDir}/quality_yield.txt",
            "~{oDir}/qualityscore.pdf",
            "~{oDir}/qualityscore.png",
            "~{oDir}/qualityscore.txt",
            "~{oDir}/sequencingArtifact.bait_bias_detail_metrics.txt",
            "~{oDir}/sequencingArtifact.bait_bias_summary_metrics.txt",
            "~{oDir}/sequencingArtifact.error_summary_metrics.txt",
            "~{oDir}/sequencingArtifact.pre_adapter_detail_metrics.txt",
            "~{oDir}/sequencingArtifact.pre_adapter_summary_metrics.txt"
        ]
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

        gpuType : select_first([runtimeAttributesOverride.acceleratorType, defaultRuntimeAttributes.acceleratorType])
        gpuCount : select_first([runtimeAttributesOverride.acceleratorCount, defaultRuntimeAttributes.acceleratorCount])

        acceleratorCount: select_first([runtimeAttributesOverride.acceleratorCount, defaultRuntimeAttributes.acceleratorCount])
        acceleratorType: select_first([runtimeAttributesOverride.acceleratorType, defaultRuntimeAttributes.acceleratorType])

        nvidiaDriverVersion : select_first([runtimeAttributesOverride.acceleratorDriverVersion, defaultRuntimeAttributes.acceleratorDriverVersion])
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : select_first([runtimeAttributesOverride.maxPreemptAttempts, defaultRuntimeAttributes.maxPreemptAttempts])
    }
}

# PARABRICKS BAMMETRICS
task bammetrics {
    input {
        File inputBAM
        File inputRefTarball

        String sampleID = ""
        String oDir = "mapping_qc/bammetrics" 


        RuntimeAttributes? runtimeAttributes

        String docker = "nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"


    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 64,
                    "cpuCount" : 16,
                    "acceleratorType" : "nvidia-tesla-a10g",
                    "acceleratorCount": 1,
                    "diskGiB" : 0,
                    "runtimeMinutes": 240,
                    "maxPreemptAttempts": 3,
                    "diskType" : "SSD",
                    "hpcQueue": "gpu",
                    "acceleratorDriverVersion": "525.60.13"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])

    # String outbase = basename(inputBAM, ".bam")
    String autoSampleID = if sampleID == "" then basename(basename(inputBAM, ".bam"), ".cram") else sampleID
    String localTarball = basename(inputRefTarball)
    String ref = basename(inputRefTarball, ".tar")

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(2.0 * size(inputBAM,  "GB")) +ceil(2.0 * size(inputRefTarball,  "GB")) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])

    


    command {
        set -euxo pipefail
        mkdir -p ~{oDir}

        ln -s ~{inputRefTarball} ${localTarball} && \
        tar xvf ~{localTarball}

        pbrun bammetrics \
        --ref ~{ref} \
        --bam ~{inputBAM} \
        --out-metrics-file ~{oDir}/~{autoSampleID}.bammetrics.txt \
        --num-threads ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])}

        rm ~{ref}*
    }

    output {
        File oMetrics= "~{oDir}/~{autoSampleID}.bammetrics.txt"
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

        gpuType : select_first([runtimeAttributesOverride.acceleratorType, defaultRuntimeAttributes.acceleratorType])
        gpuCount : select_first([runtimeAttributesOverride.acceleratorCount, defaultRuntimeAttributes.acceleratorCount])

        acceleratorCount: select_first([runtimeAttributesOverride.acceleratorCount, defaultRuntimeAttributes.acceleratorCount])
        acceleratorType: select_first([runtimeAttributesOverride.acceleratorType, defaultRuntimeAttributes.acceleratorType])

        nvidiaDriverVersion : select_first([runtimeAttributesOverride.acceleratorDriverVersion, defaultRuntimeAttributes.acceleratorDriverVersion])
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : select_first([runtimeAttributesOverride.maxPreemptAttempts, defaultRuntimeAttributes.maxPreemptAttempts])
    }
}

# SAMTOOLS COVERAGE
task samtools_coverage {
    input {
        String sampleID
        BamFile bamFile
        # File bamIndex
        
        String docker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"

        RuntimeAttributes? runtimeAttributes

    }

    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 16,
                    "cpuCount" : 6,
                    "diskGiB" : 0,
                    "runtimeMinutes": 60,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])

    String oDir = "mapping_qc/samtools_coverage"
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(2.0 * size(bamFile.bam,  "GB")) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 

    command {
        set -euxo pipefail
        mkdir -p ~{oDir}
        
        samtools coverage -o ~{oDir}/~{sampleID}.samtools_coverage.txt ~{bamFile.bam}
    }

    output {
        File coverage = "~{oDir}/~{sampleID}.samtools_coverage.txt"
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

# SOMALIER EXTRACT
task somalier_extract {
    input{
        BamFile bamFile
        # File bamFile
        # File bamIndex
        File sites
        File referenceTarball
        
        String docker ="cgrlab/somalier@sha256:3f63c21f6741e49a13dca4384840d68802eec8ead950771629bc233a055a6a1e"
        RuntimeAttributes? runtimeAttributes

    }

    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 16,
                    "cpuCount" : 8,
                    "diskGiB" : 0,
                    "runtimeMinutes": 60,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])

    String oDir = "mapping_qc/somalier/extract"
    String autoSampleID = basename(bamFile.bam, ".bam")
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(2.0 * size(bamFile.bam,  "GB")) + ceil(size(sites,  "GB"))   + ceil(2.0 *size(referenceTarball,  "GB")) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 
    String localTarball = basename(referenceTarball)
    String ref = basename(referenceTarball, ".tar")

    command {
        set -euxo pipefail
        mkdir -p ~{oDir}
        ln -s ~{referenceTarball} ${localTarball} && \
        time tar xvf ~{localTarball}

        somalier extract --sites ~{sites} -d ~{oDir} -f ~{ref}  ~{bamFile.bam}

        rm ~{ref}*
    }
    output {
        File extracted = "~{oDir}/~{autoSampleID}.somalier"
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
# SOMALIER RELATE
task somalier_relate {
    input {
        Array[File] extracted
        String project
        
        String docker ="cgrlab/somalier@sha256:3f63c21f6741e49a13dca4384840d68802eec8ead950771629bc233a055a6a1e"
        RuntimeAttributes? runtimeAttributes

    }

    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 16,
                    "cpuCount" : 2,
                    "diskGiB" : 0,
                    "runtimeMinutes": 60,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])

    String oDir = "mapping_qc/somalier/relate"
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(2.0 * size(extracted,  "GB")) + 25 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 

    command <<<
        set -euxo pipefail
        mkdir -p ~{oDir} 
        cd ~{oDir}

        somalier relate --output-prefix ~{project} ~{sep=" " extracted} 
    >>>

    output {
        File samplesTSV = "~{oDir}/~{project}.samples.tsv"
        File sampleHTML = "~{oDir}/~{project}.html"
        File ancestryPairs = "~{oDir}/~{project}.pairs.tsv"
        
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

# SOMALIER ANCESTRY
task somalier_ancestry {
    input {
        Array[File] extracted
        String outputPrefix
        Int numPCs=20
        
        File ancestryLabels
        # File labeledSamplesTarball
        String labeledSamplesDir="/data/COVID_WGS/test_qc_tools/somalier/resource/1kg-somalier"
        
        String docker ="cgrlab/somalier@sha256:3f63c21f6741e49a13dca4384840d68802eec8ead950771629bc233a055a6a1e"
        RuntimeAttributes? runtimeAttributes

    }

    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 16,
                    "cpuCount" : 2,
                    "diskGiB" : 0,
                    "runtimeMinutes": 60,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    String oDir = "mapping_qc/somalier/ancestry"
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(2.0 * size(extracted,  "GB")) + 25 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 

    command <<<
        set -euxo pipefail
        mkdir -p ~{oDir} 
        cd ~{oDir}
        somalier ancestry --n-pcs ~{numPCs} --labels ~{ancestryLabels} --output-prefix ~{outputPrefix}  ~{labeledSamplesDir}/*.somalier ++  ~{sep=" " extracted}
    >>>

    output {
        File ancestryTSV = "~{oDir}/~{outputPrefix}.somalier-ancestry.tsv"
        File ancestryHTML = "~{oDir}/~{outputPrefix}.somalier-ancestry.html"        
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

#VERIFYBAMID2
task verifybamid {
    input {
        String sampleID
        BamFile bamFile
        # File bamFile
        # File bamIndex
        File referenceTarball
        File svdTarball
        String? svdPrefix

        String docker="quay.io/biocontainers/verifybamid2:2.0.1--hbb20b25_6"

        RuntimeAttributes? runtimeAttributes

    }

    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 16,
                    "cpuCount" : 8,
                    "diskGiB" : 0,
                    "runtimeMinutes": 60,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])


    String oDir = "mapping_qc/verifybamid/~{sampleID}"
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(2.0 * size(bamFile.bam,  "GB"))  + ceil(2.0 * size(referenceTarball,  "GB")) + 50 else select_first([runtimeAttributesOverride.maxPreemptAttempts, defaultRuntimeAttributes.maxPreemptAttempts])
    String localSVDTarball = basename(svdTarball)
    String autoSVDPrefix = select_first([svdPrefix, basename(svdTarball, ".tar")] )
    String localRefTarball = basename(referenceTarball)
    String ref = basename(referenceTarball, ".tar")

    command {
        set -euxo pipefail
        ln -s ~{svdTarball} ~{localSVDTarball} && \
            time tar xvf ~{localSVDTarball}

        ln -s ~{referenceTarball} ${localRefTarball} && \
        time tar xf ~{localRefTarball}

        verifybamid2 --SVDPrefix ~{autoSVDPrefix} --Reference ~{ref} --BamFile ~{bamFile.bam} --NumThread ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])}
        mkdir -p ~{oDir}
        mv result.selfSM ~{oDir}/~{sampleID}.selfSM
        mv result.Ancestry ~{oDir}/~{sampleID}.Ancestry
        rm ~{ref}*
        rm ~{autoSVDPrefix}*
    }

    output {
        File selfSM = "~{oDir}/~{sampleID}.selfSM"
        File ancestory = "~{oDir}/~{sampleID}.Ancestry"
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

struct BamFile{
    File bam
    File bai
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