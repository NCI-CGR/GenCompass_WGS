version 1.0

workflow PremapQC {
    input {
       
        String sampleID

        File? sampleFastqLocations 
        File manifest

        File? omicsFastqMapFile
        String? awsAccountID
        
        Boolean runFastp=true
        Boolean runFastQC=true
        Boolean runFastqScreen=true
        Boolean runPostFastpQC=false

        String fastpOptions="--qualified_quality_phred 25 --n_base_limit 10 --average_qual 25 --length_required 50 --low_complexity_filter --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15"

        File? humanTarball
        File? nonhumanTarball
        Int fastqScreenSubset= 10000000
        String fastqScreenAligner="bwa"
        File? fastqScreenConfig
        File? fastqcContaminants

        String gencompassDocker=""
        String fastpDocker=""
        String fastqcDocker=""
        String fastqScreenDocker=""

        RuntimeAttributes fastpRuntimeAttributes = {
                    "memoryGiB" : 16,
                    "cpuCount" : 8
        }
        RuntimeAttributes fastqScreenRuntimeAttributes = {
                    "memoryGiB" : 16,
                    "cpuCount" : 32
        }
        RuntimeAttributes fastqcRuntimeAttributes = {
                    "memoryGiB" : 8,
                    "cpuCount" : 4
        }

        RuntimeAttributes mergeFastqRuntimeAttributes = {
                    "memoryGiB" : 8,
                    "cpuCount" : 4
        }

    }

    call prepare_paired_end_table{
        input:
            fastqLocations=sampleFastqLocations,
            manifest=manifest,
            omicsFastqMapFile=omicsFastqMapFile,
            awsAccountID=awsAccountID,
            sampleID=sampleID,
            docker=gencompassDocker,

    }

    Array[Array[String]] sampleTransformedLocations = read_tsv(prepare_paired_end_table.sampleLocations)

    scatter (lane in sampleTransformedLocations) {

        if(defined(omicsFastqMapFile)){
            call pull_omics_fastq as pullOmicsR1 {
                input:
                    omics_fastq_file=lane[2],
                    omics_fastq_map_file=select_first([omicsFastqMapFile])

            }
            call pull_omics_fastq as pullOmicsR2 {
                input:
                    omics_fastq_file=lane[3],
                    omics_fastq_map_file=select_first([omicsFastqMapFile])
            }
        }
        String laneID=lane[1]
        File r1FastqRaw= select_first([pullOmicsR1.omics_extracted_fastq_file, lane[2]])
        File r2FastqRaw = select_first([pullOmicsR2.omics_extracted_fastq_file, lane[3]])

        if(runFastp){ 
            call fastp {
                input:
                    r1Fastq=r1FastqRaw,
                    r2Fastq=r2FastqRaw,
                    sampleID=sampleID,
                    laneID=laneID,
                    fastpOptions=fastpOptions,
                    docker=fastpDocker,
                    runtimeAttributes=fastpRuntimeAttributes
            }
        }
        File r1FastpFastq = select_first([fastp.fastpR1, r1FastqRaw])
        File r2FastpFastq = select_first([fastp.fastpR2, r2FastqRaw])
        
        if(runFastQC){
            call fastqc  {
                input:
                    r1Fastq=r1FastqRaw,
                    r2Fastq=r2FastqRaw,
                    sampleID=sampleID,
                    contaminants=fastqcContaminants,
                    docker=fastqcDocker,
                    runtimeAttributes=fastqcRuntimeAttributes
             }

        }

        if(runPostFastpQC && runFastp && runFastQC){
            call fastqc as fastqc_filtered {
                input:
                    r1Fastq=r1FastpFastq,
                    r2Fastq=r2FastpFastq,
                    sampleID=sampleID,
                    contaminants=fastqcContaminants,
                    docker=fastqcDocker,
                    runtimeAttributes=fastqcRuntimeAttributes
            }
        }
    }
    if(runFastqScreen){
        File humanTarball1 = select_first([humanTarball])
        File nonhumanTarball1 = select_first([nonhumanTarball])
        call mergeFastq {
            input:
                fastq=r1FastqRaw,
                sampleID=sampleID,
                docker=gencompassDocker,
                runtimeAttributes=mergeFastqRuntimeAttributes
        }

        call fastq_screen {
            input:
                fastq = mergeFastq.mergedFastq,
                sampleID = sampleID,
                humanTarball=humanTarball1,
                nonhumanTarball=nonhumanTarball1,
                config=fastqScreenConfig,
                aligner=fastqScreenAligner,
                subset=fastqScreenSubset,
                docker=fastqScreenDocker,
                runtimeAttributes=fastqScreenRuntimeAttributes
        }

    }

    if(runPostFastpQC && runFastp &&runFastqScreen){
        File humanTarball2 = select_first([humanTarball])
        File nonhumanTarball2 = select_first([nonhumanTarball])

        call mergeFastq as mergeFastqFiltered {
        input:
            fastq=r1FastpFastq,
            sampleID=sampleID+"_filtered",
            docker=gencompassDocker,
            runtimeAttributes=mergeFastqRuntimeAttributes
        }

        call fastq_screen as fastq_screen_filtered {
            input:
                fastq = mergeFastqFiltered.mergedFastq,
                sampleID = sampleID,
                humanTarball=humanTarball2,
                nonhumanTarball=nonhumanTarball2,
                config=fastqScreenConfig,
                aligner=fastqScreenAligner,
                subset=fastqScreenSubset,
                docker=fastqScreenDocker,
                runtimeAttributes=fastqScreenRuntimeAttributes
        }

    }
    
    
    output{
        Array[File?] fastpR1 = fastp.fastpR1
        Array[File?] fastpR2 = fastp.fastpR2
        Array[File?] fastpHTML=fastp.fastpHTML
        Array[File?] fastpJSON=fastp.fastpJSON
        Array[File?] fastqcR1HTML = fastqc.r1HTML
        Array[File?] fastqcR2HTML = fastqc.r2HTML
        Array[File?] fastqcR1Zip = fastqc.r1Zip
        Array[File?] fastqcR2Zip = fastqc.r2Zip
        File? fastq_screen_HTML=fastq_screen.screenHTML
        File? fastq_screen_TXT=fastq_screen.screenTXT
        File? fastq_screen_PNG=fastq_screen.screenPNG
        File? fastq_screen_FASTQ=fastq_screen.screenFASTQ

        Array[File?] fastq_screen_filtered_R1HTML = fastqc_filtered.r1HTML
        Array[File?] fastq_screen_filtered_R2HTML = fastqc_filtered.r2HTML
        Array[File?] fastq_screen_filtered_R1Zip = fastqc_filtered.r1Zip
        Array[File?] fastq_screen_filtered_R2Zip = fastqc_filtered.r2Zip
        File? fastq_screen_filtered_HTML=fastq_screen_filtered.screenHTML
        File? fastq_screen_filtered_TXT=fastq_screen_filtered.screenTXT
        File? fastq_screen_filtered_PNG=fastq_screen_filtered.screenPNG
        File? fastq_screen_filtered_FASTQ=fastq_screen_filtered.screenFASTQ
        
    }
}


task pull_omics_fastq {
    input {
        File omics_fastq_file
        File  omics_fastq_map_file 
    }
    # Map[String, String] omics_fastq_map = read_map(omics_fastq_map_file)
    String prefix = basename(basename(omics_fastq_file, ".fastq.gz"), ".bin")
    

    command <<<
        set -eoux pipefail
        echo ~{prefix}

        python <<CODE
        import os
        import subprocess
        omics_map={}
        input_file="~{omics_fastq_file}"
        print(input_file)

        for line in open("~{omics_fastq_map_file}"):
            omics_prefx, fastq_filename = line.strip().split('\t')
            omics_map[omics_prefx] = fastq_filename
        input_prefix = os.path.basename(input_file).split('.')[0]
        output_filename = omics_map.get(input_prefix)

        print(input_prefix)
        print(omics_map)

        print(input_file)
        print(output_filename)
        subprocess.call(["cp", input_file, output_filename])
        CODE
    >>>

    runtime {
        docker: "442667281741.dkr.ecr.us-east-1.amazonaws.com/gencompass:v1.0"
        cpu : 2
        memory : "10 GiB"
    }

    output {
        File omics_extracted_fastq_file = glob("*.fastq.gz")[0]
    }
}  

task prepare_paired_end_table {
    input {
        
        File manifest
        String sampleID
        File? fastqLocations
        File? omicsFastqMapFile
        String? awsAccountID

        String docker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        Int gbRAM = 10
        
        Int runtimeMinutes = 20
        Int maxPreemptAttempts = 3
    }
    Int autoDiskGB = 10
    String oFile= "~{sampleID}.fastq_locations.tsv"

    
    command <<<
        set -euxo pipefail
        prepare_paired_end_table.py \
        ~{"--fastq_locations " + fastqLocations} \
        --manifest ~{manifest} \
        --sample_id ~{sampleID} \
        ~{"--omics_fastq_map_file " + omicsFastqMapFile } \
        ~{"--aws_account_id " + awsAccountID} \
        -odir ./
    >>>

    output {
        File sampleLocations = "~{oFile}"
    }

    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} SSD"
        cpu : 1
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        runtime_minutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
    
}

task fastp {
    input {
        File r1Fastq
        File r2Fastq
        String sampleID
        String laneID
        String fastpOptions = "--qualified_quality_phred 25 --n_base_limit 10 --average_qual 25 --length_required 50 --low_complexity_filter --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15"
        String docker = "cgrlab/fastp@sha256:a395439db5b388e0d9a97ef2c872ccd0f62af864016958cf5c42d6b411edf10b" # fastp:v0.23.2
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

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(size(r1Fastq, "GB")) + ceil(size(r2Fastq, "GB"))  + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    String oDir = "premap_qc/fastp/~{sampleID}"
    
    String r1Filename = "~{oDir}/${laneID}_R1_fastp.fastq.gz"
    String r2Filename = "~{oDir}/${laneID}_R2_fastp.fastq.gz"
    String htmlFilename = "~{oDir}/${laneID}_fastp.html"
    String jsonFilename = "~{oDir}/${laneID}_fastp.json"
    String reportTitle = "'~{sampleID}: ~{laneID}'"

    command {
        set -euxo pipefail

        mkdir -p ~{oDir}

        fastp \
        --in1 ${r1Fastq} \
        --in2 ${r2Fastq} \
        --out1 ${r1Filename} \
        --out2 ${r2Filename} \
        --report_title ${reportTitle} \
        --json ${jsonFilename} \
        --html ${htmlFilename} \
        --thread ${select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} \
        ${fastpOptions}
    }

    output {
        File fastpR1="~{r1Filename}"
        File fastpR2="~{r2Filename}"
        File fastpJSON="~{jsonFilename}"
        File fastpHTML="~{htmlFilename}"
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
        runtime_minutes  : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : select_first([runtimeAttributesOverride.maxPreemptAttempts, defaultRuntimeAttributes.maxPreemptAttempts])
    }
}

task fastqc {
    input {
        File r1Fastq
        File r2Fastq
        File? contaminants
        String sampleID

        String docker="cgrlab/qctools@sha256:d0d9514efdce0656d651a6127495bfe39a79c145714f815a338d93242eb01abf" # v2.0
        RuntimeAttributes? runtimeAttributes

    }

    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 8,
                    "cpuCount" : 4,
                    "diskGiB" : 0,
                    "runtimeMinutes": 60,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    String r1Basename = basename(basename(basename(r1Fastq, ".gz"), ".fastq"), ".fq")
    String r2Basename = basename(basename(basename(r2Fastq, ".gz"), ".fastq"), ".fq")
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(size(r1Fastq, "GiB")) + ceil(size(r2Fastq, "GiB")) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    String oDir = "premap_qc/fastqc/~{sampleID}"
    
    command {
        set -euxo pipefail
        mkdir -p ~{oDir}
        fastqc --outdir ~{oDir} \
        --threads ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} \
        ~{"--contaminants " + contaminants} ~{r1Fastq} ~{r2Fastq}
    }

    output {
        File r1HTML="~{oDir}/~{r1Basename}_fastqc.html"
        File r2HTML="~{oDir}/~{r2Basename}_fastqc.html"
        File r1Zip="~{oDir}/~{r1Basename}_fastqc.zip"
        File r2Zip="~{oDir}/~{r2Basename}_fastqc.zip"
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
        runtime_minutes  : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : select_first([runtimeAttributesOverride.maxPreemptAttempts, defaultRuntimeAttributes.maxPreemptAttempts])
    }
}

task fastq_screen {
    input {
        File fastq
        File? config
        File nonhumanTarball
        File humanTarball
        String sampleID
        Int subset=10000000
        String? aligner

        String docker="cgrlab/fastqc-screen@sha256:61203cdf4f47cea7e1e48578e9d682491bafa2e1573ae585032bde46a149ed11" # v0.15.2
        RuntimeAttributes? runtimeAttributes

    }

    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 32,
                    "cpuCount" : 16,
                    "diskGiB" : 0,
                    "runtimeMinutes": 240,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(size(fastq, "GB")) + ceil(size(nonhumanTarball, "GB")) + ceil(size(humanTarball, "GB"))+ 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    String fastqBasename = basename(basename(basename(fastq, ".gz"), ".fastq"), ".fq")
    String localNonhumanTarball = basename(nonhumanTarball)
    String localHumanTarball = basename(humanTarball)

    String oDir = "premap_qc/fastq_screen/~{sampleID}"
    
    command {
        set -euxo pipefail
        ln -s ~{nonhumanTarball} ~{localNonhumanTarball} \
            && tar xvf ~{localNonhumanTarball}
        
        mkdir fasta_human  \
            && ln -s ~{humanTarball} ~{localHumanTarball} \
            && tar xvf ~{localHumanTarball} -C fasta_human
            
        mkdir -p ~{oDir}
        fastq_screen ~{"--conf " + config} ~{"--aligner " + aligner} \
            --tag ~{"--subset " + subset} \
            --threads ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} \
            --outdir ~{oDir} ~{fastq}

        rm -r fasta_human
        rm -r fasta_nonhuman
        rm ~{localHumanTarball}
        rm ~{localNonhumanTarball}
        

        
    }

    output {
        File screenHTML="~{oDir}/~{fastqBasename}_screen.html"
        File screenTXT="~{oDir}/~{fastqBasename}_screen.txt"
        File screenPNG="~{oDir}/~{fastqBasename}_screen.png"
        File screenFASTQ="~{oDir}/~{fastqBasename}.tagged.fastq.gz"
        
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
        runtime_minutes  : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : select_first([runtimeAttributesOverride.maxPreemptAttempts, defaultRuntimeAttributes.maxPreemptAttempts])
    }
}

task mergeFastq {
    input {
        Array[File] fastq
        String sampleID

        String docker="ubuntu:22.04"
        RuntimeAttributes? runtimeAttributes

    }

    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 8,
                    "cpuCount" : 2,
                    "diskGiB" : 0,
                    "runtimeMinutes": 60,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(size(fastq, "GB")) + ceil(size(fastq, "GB")) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
        
    command {
        set -euxo pipefail
        cat ~{sep=" " fastq} > ~{sampleID}.fastq.gz
    }
    output{
        File mergedFastq = "~{sampleID}.fastq.gz"
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
        runtime_minutes  : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

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