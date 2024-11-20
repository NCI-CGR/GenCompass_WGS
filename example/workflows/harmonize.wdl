version 1.0

workflow Harmonize {
    input {
        

        File haplotypecallerVCF
        File haplotypecallerVCFIndex
    
        File deepvariantVCF
        File deepvariantVCFIndex

        File strelka2VCF
        File strelka2VCFIndex
        
        File referenceTarball
        String project = "GenCompass"
        String gencompassDocker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"

        String labelBlockSize="500M"

        RuntimeAttributes? normalizeRuntimeAttributes
        RuntimeAttributes? labelRuntimeAttributes
        RuntimeAttributes? mergeByVariantRuntimeAttributes
        RuntimeAttributes? mergeBySampleRuntimeAttributes
        
        # String called_variants
    }
    call normalizeLeftAlignSplit as normalizeHaplotypeCaller {
        input:
            caller="haplotypecaller",
            vcf=haplotypecallerVCF,
            vcfIndex=haplotypecallerVCFIndex,
            inputRefTarball=referenceTarball,
            docker=gencompassDocker,
            runtimeAttributes=normalizeRuntimeAttributes
    }

    call normalizeLeftAlignSplit as normalizeDeepVariant {
        input:
            caller="deepvariant",
            vcf=deepvariantVCF,
            vcfIndex=deepvariantVCFIndex,
            inputRefTarball=referenceTarball,
            docker=gencompassDocker,
            runtimeAttributes=normalizeRuntimeAttributes
    }

    call normalizeLeftAlignSplitStrelka as normalizeStrelka2 {
        input:
            caller="strelka2",
            vcf=strelka2VCF,
            vcfIndex=strelka2VCFIndex,
            inputRefTarball=referenceTarball,
            docker=gencompassDocker,
            runtimeAttributes=normalizeRuntimeAttributes
    }

    call labelWithCallers as labelHaplotypeCaller {
        input:
            caller="haplotypecaller",
            vcf=normalizeHaplotypeCaller.normalizedVCF,
            docker=gencompassDocker,
            runtimeAttributes=labelRuntimeAttributes,
            blockSize=labelBlockSize
    }

    call labelWithCallers as labelDeepVariant {
        input:
            caller="deepvariant",
            vcf=normalizeDeepVariant.normalizedVCF,
            docker=gencompassDocker,
            runtimeAttributes=labelRuntimeAttributes,
            blockSize=labelBlockSize
    }

    call labelWithCallers as labelStrelka2 {
        input:
            caller="strelka2",
            vcf=normalizeStrelka2.normalizedVCF,
            docker=gencompassDocker,
            runtimeAttributes=labelRuntimeAttributes,
            blockSize=labelBlockSize
    }

    Array[File] callerLabeledVCF = [labelHaplotypeCaller.labeledVCF, labelDeepVariant.labeledVCF, labelStrelka2.labeledVCF]
    Array[File] callerLabeledVCFI = [labelHaplotypeCaller.labeledVCFI, labelDeepVariant.labeledVCFI, labelStrelka2.labeledVCFI]

    call mergeByVariant {
        input:
            callerLabeledVCF=callerLabeledVCF,
            callerLabeledVCFI=callerLabeledVCFI, 
            docker=gencompassDocker,
            runtimeAttributes=mergeByVariantRuntimeAttributes
    }

    call mergeBySample {
        input:
            vcfMergedByVariant = mergeByVariant.mergedVariants,
            project=project,
            docker=gencompassDocker,
            runtimeAttributes=mergeBySampleRuntimeAttributes
    }

    output {
        File mergedVCF = mergeBySample.mergedVCF
        File mergedVCFI = mergeBySample.mergedVCFI
    }

}

struct BCFandIndex {
    File bcf
    File bcfIndex
}
task normalizeLeftAlignSplit {
    input {
        File inputRefTarball

        File vcf
        File vcfIndex

        String caller
        
        String docker = "cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 50,
                    "cpuCount" : 8,
                    "diskGiB" : 0,
                    "runtimeMinutes": 120,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(size(inputRefTarball,  "GB")) + ceil(size(vcf,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    Int nThreads = select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
    String ref = basename(inputRefTarball, ".tar")
    String oDir="~{caller}/ensemble"
    String oFile = "~{oDir}/~{caller}.normalized.vcf.gz"
    command {
        set -euxo pipefail
        time tar xvf ~{inputRefTarball}
        mkdir -p ~{oDir} 
        bcftools norm -f ~{ref} -m - --threads ~{nThreads} ~{vcf} -Oz -o ~{oFile}
        rm ~{ref}*
    }
    output {
        File normalizedVCF = "~{oFile}"
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
task normalizeLeftAlignSplitStrelka {
    input {
        File inputRefTarball

        File vcf
        File vcfIndex

        String caller
        
        String docker = "cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 50,
                    "cpuCount" : 8,
                    "diskGiB" : 0,
                    "runtimeMinutes": 120,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(size(inputRefTarball,  "GB")) + ceil(size(vcf,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    String ref = basename(inputRefTarball, ".tar")
    String oDir="~{caller}/ensemble"
    String oFile = "~{oDir}/~{caller}.normalized.vcf.gz"
    Int nThreads = select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
    command {
        set -euxo pipefail
        time tar xvf ~{inputRefTarball}
        mkdir -p ~{oDir} 
        bcftools norm -f ~{ref} -m - --threads ~{nThreads} ~{vcf} | bcftools view --min-ac 1 --threads ~{nThreads}  -Oz -o ~{oFile}
        rm ~{ref}
    }
    output {
        File normalizedVCF = "~{oFile}"
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

task labelWithCallers {
    input {
        File vcf
        String caller
        String docker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"

        # String tempDir="./"

        Map[String, String] callerLabelMap = {
            "deepvariant" : "DV", 
            "haplotypecaller": "HC",
            "strelka2" : "strelka2"
        }
        String blockSize="500M"
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 50,
                    "cpuCount" : 32,
                    "diskGiB" : 0,
                    "runtimeMinutes": 120,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    String vcfLabeledFilename = "~{caller}.labeled.vcf"
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(8.0 * size(vcf,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 
    String oDir ="~{caller}/ensemble"
    String tempDir="./"
    String callerLabel = callerLabelMap[caller]
    String vcfBasename = basename(vcf, ".gz")
    Int compression = 10
    Int nThreads = select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
    Int blocksize = ceil(size(vcf, "MB") * compression / select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount]))
    
    command<<<
        set -euxo pipefail

        mkdir -p ~{oDir} 
        mkdir work_separate && mkdir work_label

        #*********************************************
        # Separating header and variants
        #*********************************************
        echo "Separating header and variants"
        bcftools view --threads ~{nThreads} -I -h ~{vcf} -Ov -o header.vcf
        # zcat ~{vcf} | parallel --pipe --keep-order --block ~{blocksize}M -j ~{nThreads} "egrep -v '^#'" > variants.vcf
        echo "Running bgzip"
        bgzip -d ~{vcf} -c --threads ~{nThreads} > ~{vcfBasename}
        echo "Separating variants"
        time parallel --pipepart --keep-order --block ~{blockSize} --tmpdir work_separate -j ~{nThreads} -a ~{vcfBasename} "egrep -v '^#'" > variants.vcf
        rm -r work_separate
        rm ~{vcfBasename}
        
        #*********************************************
        # Label header and variant file with caller
        #*********************************************
        echo "Labeling variants with caller"
        # cat variants.vcf | parallel --pipe --keep-order --block ~{blocksize}M -j ~{nThreads} "cat > temp/{#}; prepend_labels_variant.sh temp/{#} ~{callerLabel}" > $saveDIR/~{oDir}/prepended.variants.vcf
        # time parallel --pipepart --keep-order --block ~{blockSize} --tmpdir ~{tempDir}/work -j ~{nThreads} -a ../variants.vcf "cat > {#}; prepend_labels_variant.sh {#} ~{callerLabel}" > $saveDIR/~{oDir}/prepended.variants.vcf
        time parallel --pipepart --keep-order --block ~{blockSize} --tmpdir work_label -j ~{nThreads} -a variants.vcf "cat > {#}; prepend_labels_variant.sh {#} ~{callerLabel}" > ~{oDir}/prepended.variants.vcf
        rm -r work_label
        prepend_labels_header.sh header.vcf ~{callerLabel} ~{oDir}
        rm variants.vcf
        rm header.vcf

        #*********************************************
        # Recombine labeled header and variants
        #*********************************************
        echo "Recombining header and variants"
        cd ~{oDir}
        time cat prepended.header.vcf prepended.variants.vcf > ~{vcfLabeledFilename}
        rm prepended.header.vcf 
        rm prepended.variants.vcf
        
        echo "Compressing and indexing labeled vcf"
        time bgzip --threads ~{nThreads} ~{vcfLabeledFilename}
        time tabix -p vcf ~{vcfLabeledFilename}.gz
    >>>
    
    output {
        File labeledVCF = "~{oDir}/~{vcfLabeledFilename}.gz"
        File labeledVCFI = "~{oDir}/~{vcfLabeledFilename}.gz.tbi"
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

task mergeByVariant {
    input{
        Array[File] callerLabeledVCF
        Array[File] callerLabeledVCFI
        
        String docker = "cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 50,
                    "cpuCount" : 12,
                    "diskGiB" : 0,
                    "runtimeMinutes": 240,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    String oDir="ensemble"
    String oFile = "ensemble/all_callers.vcf"

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(3.0 * size(callerLabeledVCF,  "GB"))  + ceil(size(callerLabeledVCFI,  "GB")*2) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    

    command {
        set -euxo pipefail
        mkdir -p ~{oDir}
        bcftools merge --force-samples --threads ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} -m none ~{sep=" " callerLabeledVCF} -Ov -o ~{oFile}
    }
    output{ 
        File mergedVariants = "~{oFile}"
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

task mergeBySample {
    input {
        File vcfMergedByVariant

        String docker = "cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        
        String project = ""

        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 50,
                    "cpuCount" : 32,
                    "diskGiB" : 0,
                    "runtimeMinutes": 240,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(3.0 * size(vcfMergedByVariant,  "GB")) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])

    
    String oDir = "ensemble"
    String oFileVCF = if project != "" then "~{project}.all_callers_merged_genotypes.vcf" else "all_callers_merged_genotypes.vcf"
    
    command<<<
        set -euxo pipefail
        mkdir -p ~{oDir}
        mkdir work
        parallel --pipepart --keep-order --block -1 -j ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} --tmpdir work -a ~{vcfMergedByVariant} "cat > {#};genotype_union.py {#}" >  ~{oDir}/~{oFileVCF}
        cd ~{oDir}
        bgzip --threads ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} ~{oFileVCF}
        tabix -p vcf ~{oFileVCF}.gz
    >>>

    output {
        File mergedVCF="~{oDir}/~{oFileVCF}.gz"
        File mergedVCFI="~{oDir}/~{oFileVCF}.gz.tbi"
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

