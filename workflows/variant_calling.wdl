version 1.0

workflow VariantCalling {
    input{
        String sampleID

        File sampleBAM
        File? sampleBAI
        File? sampleBQSR
        File referenceTarball
        File? callRegions

        String gencompassDocker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        String parabricksDocker="nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"
        String strelkaDocker=""

        Boolean runDeepVariant=true
        Boolean runHaplotypeCaller=true
        Boolean runStrelka2=true

        Boolean haplotypecallerUseBestPractices=false

        RuntimeAttributes deepvariantRuntimeAttributes = {
            "memoryGiB" : 64,
            "cpuCount" : 16,
            "acceleratorType" : "nvidia-tesla-a10g",
            "acceleratorCount": 1,
        }
        RuntimeAttributes haplotypecallerRuntimeAttributes = {
            "memoryGiB" : 64,
            "cpuCount" : 16,
            "acceleratorType" : "nvidia-tesla-a10g",
            "acceleratorCount": 1,
        }
        RuntimeAttributes strelkaRuntimeAttributes = {
            "memoryGiB" : 64,
            "cpuCount" : 32,
        }

        Array[String] deepVariantAdditionalArgs=["--num-streams-per-gpu 4",  "--run-partition",  "--gpu-num-per-partition 1"]
        Array[String] haplotypecallerAdditionalArgs=[]
    }
    if(!defined(sampleBAI)){
        call generateBamIndex {
            input:
                bamInput=sampleBAM,
                docker=gencompassDocker
        }
    }

    BamFile bamFile={
            "bam" : select_first([generateBamIndex.bam, sampleBAM]),
            "bai": select_first([generateBamIndex.bai, sampleBAI])
        }


    if(runDeepVariant){

        call deepvariant_parabricks4 as deepvariant {
            input:
                sampleID=sampleID,
                bamFile=bamFile,
                inputRefTarball=referenceTarball,
                intervalFile=callRegions,
                docker=parabricksDocker,
                runtimeAttributes=deepvariantRuntimeAttributes,
                additionalArgs=deepVariantAdditionalArgs
        }

        call compressVCF as compressDeepVariant {
            input:
                vcf = deepvariant.deepvariantVCF,
                oDir="deepvariant/~{sampleID}",
                docker=gencompassDocker
        }
    }
    
    if(runHaplotypeCaller){
        call haplotypecaller_parabricks4 as haplotypecaller {
            input:
                sampleID=sampleID,
                bamFile=bamFile,
                inputRecal=sampleBQSR,
                intervalFile=callRegions,
                inputRefTarball=referenceTarball,
                useBestPractices=haplotypecallerUseBestPractices,
                oDir="haplotypecaller/~{sampleID}",
                docker=parabricksDocker,
                additionalArgs=haplotypecallerAdditionalArgs,
                runtimeAttributes=haplotypecallerRuntimeAttributes

        }
        call compressVCF as compressHaplotypeCaller {
            input:
                vcf = haplotypecaller.haplotypecallerVCF,
                oDir="haplotypecaller/~{sampleID}",
                docker=gencompassDocker
        }
    }
    
    if(runStrelka2){
        call strelka {
            input:
                sampleID=sampleID,
                bamFile=bamFile,
                callRegions=callRegions,
                inputRefTarball=referenceTarball,
                oDir="strelka2/~{sampleID}",
                docker=strelkaDocker,
                runtimeAttributes=strelkaRuntimeAttributes
        }
    }
    
        


    output{
        
        File? deepvariantVCF=compressDeepVariant.vcfGZ
        File? deepvariantTBI=compressDeepVariant.vcfTBI

        File? haplotypecallerVCF = compressHaplotypeCaller.vcfGZ
        File? haplotypecallerTBI = compressHaplotypeCaller.vcfTBI

        File? strelkaVCF=strelka.strelkaVCF
        File? strelkaTBI = strelka.strelkaTBI
        File? strelkaGVCF=strelka.strelkaGVCF
        File? strelkaGTBI = strelka.strelkaGTBI
    }
}

task generateBamIndex {
    input {
        File bamInput
        String docker
        Int runtimeMinutes=60
        Int nThreads=8
        Int gbRAM=50
        Int diskGB=0


    }
    String localBAM = basename(bamInput)
        Int autoDiskGB = if diskGB < 1 then ceil(3.0 * size(bamInput,  "GB")) +  50 else diskGB

    command <<<
        set -euxo pipefail
        cp ~{bamInput} ./
        set -euox pipefail
        samtools index -b -@ ~{nThreads} ~{localBAM}
    >>>

    output {
        File bam=localBAM
        File bai="~{localBAM}.bai"

    }
    runtime {
        docker: "~{docker}"
        runtime_minutes : runtimeMinutes
        cpu : "~{nThreads}"
        memory : "~{gbRAM} GiB"
        diskGB : autoDiskGB
        queue : "norm"

        disks : "local-disk ~{autoDiskGB} SSD"

        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : runtimeMinutes
        preemptible : 3
    }
}
task deepvariant_parabricks4 {
    input {
        String sampleID
        BamFile bamFile
        File inputRefTarball
        File? intervalFile
        Boolean gvcfMode = true
        String docker = "nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"
        
        Array[String] additionalArgs=["--num-streams-per-gpu 4",  "--run-partition",  "--gpu-num-per-partition 1"]
        RuntimeAttributes? runtimeAttributes 
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 64,
                    "cpuCount" : 16,
                    "acceleratorType" : "nvidia-tesla-a10g",
                    "acceleratorCount": 1,
                    "diskGiB" : 0,
                    "acceleratorDriverVersion": "460.73.01",
                    "runtimeMinutes": 240,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "gpu",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    # String outbase = basename(inputBAM, ".bam")
    String localTarball = basename(inputRefTarball)
    String ref = basename(inputRefTarball, ".tar")
    String outVCF = sampleID + ".deepvariant" + (if gvcfMode then '.g' else '') + ".vcf"
    String localBAM = sampleID + ".bam"
    
    String oDir="deepvariant/~{sampleID}"
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(3.0 * size(bamFile.bam,  "GB")) + ceil(3.0 * size(inputRefTarball,  "GB")) + ceil(size(bamFile.bai,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 

    # Int mb_ram = gbRAM * 1024

    command <<<
        set -euxo pipefail
        mkdir -p ~{oDir} && \
        cd ~{oDir}
        # cp required for omics
        cp ~{bamFile.bam} ~{localBAM}
        cp ~{bamFile.bai} ~{localBAM}.bai

        ln -s ~{inputRefTarball} ~{localTarball} && \
        tar xf ~{localTarball}
        
        pbrun deepvariant \
            --num-gpu ~{select_first([runtimeAttributesOverride.acceleratorCount, defaultRuntimeAttributes.acceleratorCount])} \
            --ref ~{ref} \
            ~{if gvcfMode then "--gvcf " else ""} \
            ~{"--interval-file " + intervalFile} \
            --in-bam ~{localBAM} \
            --out-variants ~{outVCF} \
            ~{sep = " " additionalArgs}
        
        rm ~{localBAM}
        rm ~{localBAM}.bai
        rm ~{ref}*
        
        >>>

    output {
        File deepvariantVCF = "~{oDir}/~{outVCF}"
        # File deepvariantTBI = "~{oDir}/~{outVCF}.gz.tbi"
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

task haplotypecaller_parabricks4{
    input {
        BamFile bamFile
        File? inputRecal
        String sampleID
        File inputRefTarball
        File? intervalFile
        String docker = "nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1"

        String oDir="haplotypecaller"

        Boolean gvcfMode = true
        Boolean useBestPractices = false
        String? haplotypecallerPassthroughOptions
        String annotationArgs = ""
       
        Array[String] additionalArgs=[]
        RuntimeAttributes? runtimeAttributes 
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 64,
                    "cpuCount" : 16,
                    "acceleratorType" : "nvidia-tesla-a10g",
                    "acceleratorCount": 1,
                    "diskGiB" : 0,
                    "acceleratorDriverVersion": "460.73.01",
                    "runtimeMinutes": 240,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "gpu",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    # String outbase = basename(inputBAM, ".bam")
    String localTarball = basename(inputRefTarball)
    String ref = basename(inputRefTarball, ".tar")
    String localBAM = sampleID + ".bam"

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(3.0 * size(bamFile.bam,  "GB")) + ceil(3.0 * size(inputRefTarball,  "GB")) + ceil(size(bamFile.bai,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 

    String outVCF = sampleID + ".haplotypecaller" + (if gvcfMode then '.g' else '') + ".vcf"

    String quantization_band_stub = if useBestPractices then " -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 " else ""
    String quantization_qual_stub = if useBestPractices then " --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30" else ""
    String annotation_stub_base = if useBestPractices then "-G StandardAnnotation -G StandardHCAnnotation" else annotationArgs
    String annotation_stub = if useBestPractices && gvcfMode then annotation_stub_base + " -G AS_StandardAnnotation " else annotation_stub_base


    command <<<
        set -euxo pipefail
        mkdir -p ~{oDir}
        cd ~{oDir}
        cp ~{bamFile.bam} ~{localBAM}
        cp ~{bamFile.bai} ~{localBAM}.bai

        ln -s ~{inputRefTarball} ~{localTarball} && \
        tar xvf ~{localTarball}

        pbrun haplotypecaller \
        --in-bam ~{localBAM} \
        --ref ~{ref} \
        --out-variants ~{outVCF} \
        ~{"--in-recal-file " + inputRecal} \
        ~{if gvcfMode then "--gvcf " else ""} \
        ~{"--interval-file " + intervalFile} \
        ~{"--haplotypecaller-options " + '"' + haplotypecallerPassthroughOptions + '"'} \
        ~{annotation_stub} \
        ~{quantization_band_stub} \
        ~{quantization_qual_stub} \
        ~{sep = " " additionalArgs} >> parabricks_log.txt

        rm ~{ref}*
        rm ~{localBAM}
        rm ~{localBAM}.bai

    >>>

    output {
        File haplotypecallerVCF = "~{oDir}/~{outVCF}"
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
task strelka {
    input {
        String sampleID
        BamFile bamFile
        File inputRefTarball
        File? callRegions
        String docker="quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"

        String oDir="strelka2"


        String mode="local"
        RuntimeAttributes? runtimeAttributes
    }
    
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 64,
                    "cpuCount" : 32,
                    "diskGiB" : 0,
                    "runtimeMinutes": 240,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    # String outbase = basename(inputBAM, ".bam")
    String ref = basename(inputRefTarball, ".tar")
    String localBAM = sampleID + ".bam"

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(3.0 * size(bamFile.bam,  "GB")) + ceil(3.0 * size(inputRefTarball,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 

    # Int mb_ram = gbRAM * 1024

    command <<<
        set -euox pipefail
        
        mkdir -p ~{oDir}

        cp ~{bamFile.bam} ~{localBAM}
        cp ~{bamFile.bai} ~{localBAM}.bai
        # cp ~{bamFile.bam} ./
        # cp ~{bamFile.bai} ./
        time tar xvf ~{inputRefTarball}
        
        ~{"bgzip -c " + callRegions + " > callregions.bed.gz; tabix callregions.bed.gz"}

        configureStrelkaGermlineWorkflow.py \
        --bam ~{localBAM} \
        --referenceFasta ~{ref} \
        ~{if defined(callRegions) then "--callRegions callregions.bed.gz" else ""} \
        --runDir ./

        ./runWorkflow.py \
        -m ~{mode} \
        -j ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} \
        -g ~{select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB])}

        cp results/variants/variants.vcf.gz ~{oDir}/~{sampleID}.strelka.variants.vcf.gz
        cp results/variants/variants.vcf.gz.tbi ~{oDir}/~{sampleID}.strelka.variants.vcf.gz.tbi 
        cp results/variants/genome.vcf.gz ~{oDir}/~{sampleID}.strelka.genome.vcf.gz 
        cp results/variants/genome.vcf.gz.tbi ~{oDir}/~{sampleID}.strelka.genome.vcf.gz.tbi 

        rm ~{ref}*
        rm ~{localBAM}
        rm ~{localBAM}.bai

    >>>

    output {
        File strelkaVCF = "~{oDir}/~{sampleID}.strelka.variants.vcf.gz"
        File strelkaTBI = "~{oDir}/~{sampleID}.strelka.variants.vcf.gz.tbi"
        File strelkaGVCF= "~{oDir}/~{sampleID}.strelka.genome.vcf.gz"
        File strelkaGTBI="~{oDir}/~{sampleID}.strelka.genome.vcf.gz.tbi"
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
task compressVCF {
    input {
        File vcf
        String oDir
        String docker = "cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"

        Int nThreads = 8
        Int gbRAM = 5
    }
    Int autoDiskGB = ceil(size(vcf, "GiB")*2.5)
    String vcfBasename = basename(vcf)

    command {
        set -euxo pipefail
        mkdir -p ~{oDir} && cd ~{oDir}
        ln -s ~{vcf} ./

        bgzip --threads ~{nThreads} ~{vcfBasename}
        tabix ~{vcfBasename}.gz
    }

    runtime {
        docker : "~{docker}"
        disks : "local-disk ~{autoDiskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : 60
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : 3
    }

    output {
        File vcfGZ = "~{oDir}/~{vcfBasename}.gz"
        File vcfTBI = "~{oDir}/~{vcfBasename}.gz.tbi"
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
