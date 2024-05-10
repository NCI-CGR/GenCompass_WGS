version 1.0

struct BamFile{
    File bam
    File bai
}

task extract {
    input{
        BamFile bamFile
        # File bamFile
        # File bamIndex
        File sites
        File referenceTarball
        
        String docker ="cgrlab/somalier@sha256:3f63c21f6741e49a13dca4384840d68802eec8ead950771629bc233a055a6a1e"
        Int nThreads=8
        String hpcQueue="norm"
        Int gbRAM=2
        Int diskGB=0
        Int runtimeMinutes=60
        Int maxPreemptAttempts=3
    }

    String oDir = "mapping_qc/somalier/extract"
    String autoSampleID = basename(bamFile.bam, ".bam")
    Int auto_diskGB = if diskGB < 1 then ceil(2.0 * size(bamFile.bam,  "GiB"))  + ceil(size(referenceTarball,  "GiB")) + 50 else diskGB
    String localTarball = basename(referenceTarball)
    String ref = basename(referenceTarball, ".tar")

    command {
        set -ex
        mkdir -p ~{oDir}
        mv ~{referenceTarball} ${localTarball} && \
        time tar xvf ~{localTarball}

        somalier extract --sites ~{sites} -d ~{oDir} -f ~{ref}  ~{bamFile.bam}

        rm ~{ref}*
    }
    output {
        File extracted = "~{oDir}/~{autoSampleID}.somalier"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM}  GiB"
        hpcMemory : gbRAM
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

task relate {
    input {
        Array[File] extracted
        String project
        
        String docker ="cgrlab/somalier@sha256:3f63c21f6741e49a13dca4384840d68802eec8ead950771629bc233a055a6a1e"
        Int nThreads=8
        String hpcQueue="norm"
        Int gbRAM=2
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
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM}  GiB"
        hpcMemory : gbRAM
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

task ancestry {
    input {
        Array[File] extracted
        Int numPCs=20
        String project
        
        File ancestryLabels
        # File labeledSamplesTarball
        String labeledSamplesDir="/data/COVID_WGS/test_qc_tools/somalier/resource/1kg-somalier"
        
        String docker ="cgrlab/somalier@sha256:3f63c21f6741e49a13dca4384840d68802eec8ead950771629bc233a055a6a1e"
        Int nThreads=8
        String hpcQueue="norm"
        Int gbRAM=2
        Int diskGB=0
        Int runtimeMinutes=60
        Int maxPreemptAttempts=3

    }
    String oDir = "mapping_qc/somalier/ancestry"
    Int auto_diskGB = if diskGB < 1 then 500 else diskGB
    # String localTarball = basename(labeledSamplesTarball)
    # String labeledSamplesDir = basename(labeledSamplesTarball, ".tar")
    # mv ~{labeledSamplesTarball} ${localTarball}
    # tar xvf ~{localTarball}
    command <<<
        mkdir -p ~{oDir} 
        cd ~{oDir}
        somalier ancestry --n-pcs ~{numPCs} --labels ~{ancestryLabels} --output-prefix ~{project}  ~{labeledSamplesDir}/*.somalier ++  ~{sep=" " extracted}
    >>>

    output {
        File ancestryTSV = "~{oDir}/~{project}.somalier-ancestry.tsv"
        File ancestryHTML = "~{oDir}/~{project}.somalier-ancestry.html"        
    }
    runtime {
        docker : docker
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM}  GiB"
        hpcMemory : gbRAM
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

workflow Somalier {
    input {
        String somalierDocker ="cgrlab/somalier@sha256:3f63c21f6741e49a13dca4384840d68802eec8ead950771629bc233a055a6a1e"
        File somalierExtractFOF
        Boolean runRelate=true
        Boolean runAncestry=false
        String project
        File? ancestryLabels
    }
    Array[File] somalierExtract = read_lines(somalierExtractFOF)
    if(runRelate){
        call relate{
            input:
                docker=somalierDocker,
                extracted=somalierExtract,
                project=project
        }
    }
    if(runAncestry) {
        File ancestryLabelsX = select_first([ancestryLabels])
        call ancestry{
            input:
            docker=somalierDocker,
            extracted=somalierExtract,
            ancestryLabels=ancestryLabelsX,
            project=project

        }
    }
}