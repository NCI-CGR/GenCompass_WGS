version 1.0

workflow Harmonize {
    input {

        File intervalSplitFile
        

        File haplotypecallerVCF
        File haplotypecallerVCFIndex
    
        File deepvariantVCF
        File deepvariantVCFIndex

        File strelka2VCF
        File strelka2VCFIndex
        
        # File referenceTarball
        Reference reference
        String project = "GenCompass"
        String gencompassDocker
        String labelBlockSize="-2"
        String mergeBySampleBlockSize="-2"

        RuntimeAttributes? normalizeRuntimeAttributes
        RuntimeAttributes? mergeByVariantRuntimeAttributes
        RuntimeAttributes? mergeBySampleRuntimeAttributes
        RuntimeAttributes? concatSplitIntervalsRuntimeAttributes

        # String called_variants
    }
    Array[Array[String]] intervalSplitRegions = read_tsv(intervalSplitFile)

    scatter(intervalSplit in intervalSplitRegions) {
        String intervalSplitName = intervalSplit[0] + "." + intervalSplit[1] + "." + intervalSplit[2]
        String intervalSplitRegion = intervalSplit[0] + ":" + intervalSplit[1] + "-" + intervalSplit[2]

        call splitNormalizeLabelHaplotypeCallerDeepVariant as normalizeHaplotypeCaller {
            input:
                project="~{project}.~{intervalSplitName}",
                caller="haplotypecaller",
                interval=intervalSplitRegion,
                vcf=haplotypecallerVCF,
                vcfIndex=haplotypecallerVCFIndex,
                reference=reference,
                blockSize=labelBlockSize,
                docker=gencompassDocker,
                runtimeAttributes=normalizeRuntimeAttributes
        }

        call splitNormalizeLabelHaplotypeCallerDeepVariant as normalizeDeepVariant {
            input:
                project="~{project}.~{intervalSplitName}",
                caller="deepvariant",
                interval=intervalSplitRegion,
                vcf=deepvariantVCF,
                vcfIndex=deepvariantVCFIndex,
                reference=reference,
                blockSize=labelBlockSize,
                docker=gencompassDocker,
                runtimeAttributes=normalizeRuntimeAttributes
        }

        call splitNormalizeLabelStrelka as normalizeStrelka2 {
            input:
                project="~{project}.~{intervalSplitName}",
                caller="strelka2",
                interval=intervalSplitRegion,
                vcf=strelka2VCF,
                vcfIndex=strelka2VCFIndex,
                blockSize=labelBlockSize,
                reference=reference,
                docker=gencompassDocker,
                runtimeAttributes=normalizeRuntimeAttributes
        }
        
        Array[File] callerLabeledVCF = [normalizeHaplotypeCaller.labeledVCF, normalizeDeepVariant.labeledVCF, normalizeStrelka2.labeledVCF]
        Array[File] callerLabeledVCFI = [normalizeHaplotypeCaller.labeledVCFI, normalizeDeepVariant.labeledVCFI, normalizeStrelka2.labeledVCFI]

        call mergeByVariant {
            input:
                project="~{project}.~{intervalSplitName}",
                callerLabeledVCF=callerLabeledVCF,
                callerLabeledVCFI=callerLabeledVCFI, 
                docker=gencompassDocker,
                runtimeAttributes=mergeByVariantRuntimeAttributes
        }

        call mergeBySample {
            input:
                vcfMergedByVariant = mergeByVariant.mergedVariants,
                blockSize=mergeBySampleBlockSize,
                project="~{project}.~{intervalSplitName}",
                docker=gencompassDocker,
                runtimeAttributes=mergeBySampleRuntimeAttributes
        }
    }
    call concatSplitIntervals {
        input:
            harmonizedSplitVCF=mergeBySample.mergedVCF,
            harmonizedSplitVCFI=mergeBySample.mergedVCFI,
            project=project,
            docker=gencompassDocker,
            runtimeAttributes=concatSplitIntervalsRuntimeAttributes


    }

    output {
        File harmonizedVCF = concatSplitIntervals.harmonizedMergedVCF
        File harmonizedVCFI = concatSplitIntervals.harmonizedMergedVCFI
        
    }
}


struct BCFandIndex {
    File bcf
    File bcfIndex
}


task splitVariantByInterval {
    input {
        File vcf
        File vcfIndex
        String interval
        String name
        String docker
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 25,
                    "cpuCount" : 16,
                    "diskGiB" : 0,
                    "runtimeMinutes": 180,
                    "hpcQueue": "norm",
                    "diskType" : "HDD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    
    Int nThreads = select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
    # String oDir = "splitVariants"
    
    command <<<
        set -euxo pipefail
        bcftools filter --regions ~{interval} --write-index -Ob -o ~{name}.bcf.gz ~{vcf}
    >>>

    output {
        File filteredVCF="~{name}.bcf.gz"
        File filteredVCFI="~{name}.vcf.gz.csi"
    }
    runtime {
        docker : docker
        cpu : select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
        memory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) + " GiB"
        
        hpcMemory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB])
        memory_mb : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) * 1024
        hpcQueue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        queue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        hpcRuntimeMinutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])
        runtime_minutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

    }
}

task splitNormalizeLabelStrelka {
    input {
        File vcf
        File vcfIndex
        String interval
        String project
        String caller
        String docker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        Reference reference

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
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    String vcfLabeledFilename = "~{project}.~{caller}.labeled.vcf"
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(8.0 * size(vcf,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 
    String oDir ="~{caller}/ensemble"
    String callerLabel = callerLabelMap[caller]
    Int nThreads = select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
    
    command<<<
        set -euxo pipefail
        mkdir -p ~{oDir} 


        #*********************************************
        # Update header
        #*********************************************
        bcftools view --header-only --no-update ~{vcf} -o header.txt
        sed -i 's\^##FORMAT=<ID=AD,Number=.\##FORMAT=<ID=AD,Number=R\g' header.txt

        #*********************************************
        # Normalize
        #   - reheader
        #   - filter based on call regions
        #   - normalize
        #   - filter min allele count 1
        #*********************************************
        echo "Normalizing"
        bcftools reheader --header header.txt ~{vcf} --threads ~{nThreads} -o reheaded.vcf
        bcftools index reheaded.vcf  # optional but recommended if compressed


        bcftools filter --threads ~{nThreads} --regions ~{interval} -Ou reheaded.vcf | \
        bcftools norm -f ~{reference.fasta} -m - --threads ~{nThreads} -Ou | \
        bcftools view --min-ac 1 --threads ~{nThreads}  -Ov -o normalized.vcf

        rm reheaded.vcf

        #*********************************************
        # Label vcf with caller
        #*********************************************
        echo "Labeling variants with caller"
        mkdir work_label
        time parallel --pipepart --keep-order --block ~{blockSize} --tmpdir work_label -j ~{nThreads} -a normalized.vcf "prepend_labels_stdin.sh ~{callerLabel}" > ~{oDir}/~{vcfLabeledFilename}
        # rm -r work_label
        rm normalized.vcf
        cd ~{oDir}


        #*********************************************
        # Compress
        #*********************************************
        echo "Compressing and indexing labeled vcf"
        time bgzip --threads ~{nThreads} ~{vcfLabeledFilename}
        time bcftools index --threads ~{nThreads}  ~{vcfLabeledFilename}.gz
    >>>
    
    output {
        File labeledVCF = "~{oDir}/~{vcfLabeledFilename}.gz"
        File labeledVCFI = "~{oDir}/~{vcfLabeledFilename}.gz.csi"
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
    }

}

task splitNormalizeLabelHaplotypeCallerDeepVariant {
    input {
        File vcf
        File vcfIndex
        String interval
        Reference reference
        String project
        String caller
        String docker

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
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    String vcfLabeledFilename = "~{project}.~{caller}.labeled.vcf"
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])  < 1 then ceil(8.0 * size(vcf,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 
    String oDir ="~{caller}/ensemble"
    String callerLabel = callerLabelMap[caller]
    Int nThreads = select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
    
    command<<<
        set -euxo pipefail
        mkdir -p ~{oDir} 
        

        #*********************************************
        # Normalize
        #*********************************************
        echo "Normalizing"
        time bcftools filter --threads ~{nThreads} --regions ~{interval} -Ou ~{vcf} | \
            bcftools norm --threads ~{nThreads} -f ~{reference.fasta} -m - -Ov -o normalized.vcf

        #*********************************************
        # Label header and variant file with caller
        #*********************************************
        echo "Labeling variants with caller"
        mkdir work_label
        time parallel --pipepart --keep-order --block ~{blockSize} --tmpdir work_label -j ~{nThreads} -a normalized.vcf "prepend_labels_stdin.sh ~{callerLabel}" > ~{oDir}/~{vcfLabeledFilename}
        rm normalized.vcf
        cd ~{oDir}
        
        # echo "Normalizing and labelling vcf"
        # time bcftools filter --regions ~{interval} -Ou ~{vcf} | \
        #     bcftools norm -f ~{reference.fasta} -m - --threads ~{nThreads} -Ov | \
        #     prepend_labels_stdin.sh ~{callerLabel} | \
        #     bcftools view -Oz --write-index --threads ~{nThreads} -o ~{vcfLabeledFilename}.gz

        #*********************************************
        # Compress normalized + labeled VCF
        #*********************************************
        echo "Compressing and indexing labeled vcf"
        time bgzip --threads ~{nThreads} ~{vcfLabeledFilename}
        time bcftools index --threads ~{nThreads}  ~{vcfLabeledFilename}.gz
    >>>
    
    output {
        File labeledVCF = "~{oDir}/~{vcfLabeledFilename}.gz"
        File labeledVCFI = "~{oDir}/~{vcfLabeledFilename}.gz.csi"
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
    }

}

task mergeByVariant {
    input{
        Array[File] callerLabeledVCF
        Array[File] callerLabeledVCFI
        String project
        
        String docker
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 50,
                    "cpuCount" : 12,
                    "diskGiB" : 0,
                    "runtimeMinutes": 240,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    String oDir="ensemble"
    String oFile = "ensemble/~{project}.all_callers.vcf"

    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(3.0 * size(callerLabeledVCF,  "GB"))  + ceil(size(callerLabeledVCFI,  "GB")*2) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    Int nThreads = select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])

    command {
        set -euxo pipefail
        mkdir -p ~{oDir}
        bcftools merge --force-samples --threads ~{nThreads} -m none ~{sep=" " callerLabeledVCF} -Ov -o ~{oFile}
    }
    output{ 
        File mergedVariants = "~{oFile}"
    }

    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} ~{select_first([runtimeAttributesOverride.diskType, defaultRuntimeAttributes.diskType])}"
        cpu : nThreads
        memory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) + " GiB"
        
        hpcMemory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB])
        memory_mb : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) * 1024
        hpcQueue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        queue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        hpcRuntimeMinutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])
        runtime_minutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

    }
    

}

task mergeBySample {
    input {
        File vcfMergedByVariant

        String docker
        
        String project = ""
        String blockSize="-1"

        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 50,
                    "cpuCount" : 32,
                    "diskGiB" : 0,
                    "runtimeMinutes": 240,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(3.0 * size(vcfMergedByVariant,  "GB")) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])

    
    String oDir = "ensemble"
    String oFileVCF = if project != "" then "~{project}.all_callers_merged_genotypes.vcf" else "all_callers_merged_genotypes.vcf"
    Int nThreads = select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
    
    command<<<
        set -euxo pipefail
        mkdir -p ~{oDir}
        mkdir -p work/work
        cd work
        parallel --pipepart --keep-order --block ~{blockSize} -j ~{nThreads} --tmpdir work -a ~{vcfMergedByVariant} "cat > {#};genotype_union.py {#}" >  ../~{oDir}/~{oFileVCF}
        cd ../
        rm -r work
        cd ~{oDir}
        bgzip --threads ~{nThreads} ~{oFileVCF}
        bcftools index --threads ~{nThreads}  ~{oFileVCF}.gz
    >>>

    output {
        File mergedVCF="~{oDir}/~{oFileVCF}.gz"
        File mergedVCFI="~{oDir}/~{oFileVCF}.gz.csi"
    }

    runtime {
        docker : docker
        disks : "local-disk ~{autoDiskGB} ~{select_first([runtimeAttributesOverride.diskType, defaultRuntimeAttributes.diskType])}"
        cpu : nThreads
        memory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) + " GiB"
        
        hpcMemory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB])
        memory_mb : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) * 1024
        hpcQueue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        queue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        hpcRuntimeMinutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])
        runtime_minutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

    }

}

task concatSplitIntervals {
    input{
        Array[File] harmonizedSplitVCF
        Array[File] harmonizedSplitVCFI
        String project
        String docker
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                "memoryGiB" : 50,
                "cpuCount" : 8,
                "runtimeMinutes": 120,
                "hpcQueue": "norm",
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int nThreads = select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])
    String oDir = "ensemble"
    String oFileVCF = if project != "" then "~{project}.all_callers_merged_genotypes.vcf.gz" else "all_callers_merged_genotypes.vcf.gz"
    command <<<
        set -euxo pipefail
        mkdir -p ~{oDir} 
        bcftools concat --threads ~{nThreads} --write-index -Oz -o ~{oDir}/~{oFileVCF} ~{sep=' ' harmonizedSplitVCF}
    >>>
    output{
        File harmonizedMergedVCF="~{oDir}/~{oFileVCF}"
        File harmonizedMergedVCFI="~{oDir}/~{oFileVCF}.csi"
    }
    runtime {
        docker : docker
        cpu : nThreads
        memory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) + " GiB"
        
        hpcMemory : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB])
        memory_mb : select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB]) * 1024
        hpcQueue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        queue : select_first([runtimeAttributesOverride.hpcQueue, defaultRuntimeAttributes.hpcQueue])
        hpcRuntimeMinutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])
        runtime_minutes : select_first([runtimeAttributesOverride.runtimeMinutes, defaultRuntimeAttributes.runtimeMinutes])

    }
}

struct Reference {
    File fasta
    Array[File] index
}


struct RuntimeAttributes {
    Int? memoryGiB
    Int? cpuCount
    Int? diskGiB
    Int? bootDiskGiB
    
    Int? runtimeMinutes
    String? hpcQueue
    String? diskType
}

