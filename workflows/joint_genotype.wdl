version 1.0

workflow JointGenotype {
    input{
        File? variantFOF 
        Array[File] calledVariants =[]
        String caller
        File intervalBED
        Int nucleotidesPerBin=85000000
        String gencompassDocker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"

        String glnexusDocker

        Boolean splitIntervals=true
        Boolean skipCombine=true
        Boolean keepIndividualBCFs=true

        String glnexusConfig

        RuntimeAttributes? binIntervalsRuntimeAttributes
        RuntimeAttributes? splitIntervalsRuntimeAttributes
        RuntimeAttributes? glnexusRuntimeAttributes
        RuntimeAttributes? concatRuntimeAttributes
        RuntimeAttributes? indexRuntimeAttributes
        RuntimeAttributes? compressRuntimeAttributes


    }
    
    String glnexusBaseConfig = basename(glnexusConfig, ".yml")
    if(glnexusBaseConfig == glnexusConfig){
        String glnexusConfigOption = glnexusConfig
    }
    if(glnexusBaseConfig != glnexusConfig){
        File glnexusConfigFile=glnexusConfig
    }

    Array[File] variantFiles = read_lines(select_first([variantFOF, write_lines(calledVariants)])) 

    if(splitIntervals){
        call binIntervals{
            input:
                bedFile=intervalBED,
                docker=gencompassDocker,
                nucleotidesPerBin=nucleotidesPerBin,
                runtimeAttributes=binIntervalsRuntimeAttributes
        }

        scatter(variantPath in variantFiles) {
            File variantVCF=variantPath
            File variantIndex=variantPath+".tbi"
            Variant variant = {
                "variant": variantVCF, 
                "variantIndex": variantIndex
            }
            call splitVariantByIntervalBins {
                input:
                    variant=variant,
                    intervalBins=binIntervals.binnedIntervals,
                    docker=gencompassDocker,
                    runtimeAttributes=splitIntervalsRuntimeAttributes
            }
        }

        Array[Array[File]] transposedVariantIntervals = transpose(splitVariantByIntervalBins.variantInterval)

        scatter(binnedIntervalPair in zip(binIntervals.binnedIntervals, transposedVariantIntervals)){
            String intervalName=basename(binnedIntervalPair.left, ".bed")
            call glnexus {
                input:
                    variants=binnedIntervalPair.right,
                    intervalName=intervalName,
                    caller=caller,
                    glnexusConfigOption=glnexusConfigOption,
                    glnexusConfigFile=glnexusConfigFile,
                    docker=glnexusDocker,
                    runtimeAttributes=glnexusRuntimeAttributes
                }
                if(keepIndividualBCFs){
                    call compressBCF as compressIndividualBCF {
                    input:
                        oDir="joint_genotype/~{caller}",
                        bcf=glnexus.glnexusBCF,
                        docker=gencompassDocker,
                        runtimeAttributes=compressRuntimeAttributes
                    }
                }
                if(!skipCombine){
                    call indexBCF as indexGlnexus{
                        input:
                            bcf=glnexus.glnexusBCF,
                            oDir="joint_genotype/~{caller}",
                            docker=gencompassDocker
                    }
                }
                File bcfFile = select_first([indexGlnexus.bcfOriginal, glnexus.glnexusBCF])
                File bcfIndexFile = select_first([indexGlnexus.bcfIndex, glnexus.glnexusBCF])


            }

            BcfAndIndexArray bcfsAndIndices = {
                "bcf": bcfFile,
                "bcfIndex": bcfIndexFile
            }
            if(!skipCombine){
                call concatVariantIntervals {
                input:
                    bcfsAndIndices = bcfsAndIndices,
                    caller=caller,
                    docker=gencompassDocker,
                    runtimeAttributes=concatRuntimeAttributes
            }
        }
    }
    if(!splitIntervals) {
        String fullIntervalName=basename(intervalBED, ".bed")
        call glnexus as glnexus_single {

            input:
                variants=variantFiles,
                intervalName=fullIntervalName,
                intervalBED=intervalBED,
                caller=caller,
                glnexusConfigOption=glnexusConfigOption,
                glnexusConfigFile=glnexusConfigFile,
                docker=glnexusDocker,
                runtimeAttributes=indexRuntimeAttributes
        }
    }

    if(!skipCombine || !splitIntervals){
        File glnexusBCF = select_first([concatVariantIntervals.concatVariants, glnexus_single.glnexusBCF])

        call compressBCF {
        input:
            oDir="joint_genotype/~{caller}",
            bcf=glnexusBCF,
            docker=gencompassDocker,
            runtimeAttributes=compressRuntimeAttributes
        }

    }
    
    output{
        Array[File?]? individualBCF = compressIndividualBCF.bcfCompressed
        Array[File?]? individualBCFIndex = compressIndividualBCF.bcfIndex
        # BcfAndIndexArray? bcfAndIndexArray = bcfsAndIndices
        File? concatVariants = compressBCF.bcfCompressed
        File? concatVariantsIndex = compressBCF.bcfIndex
        File? binnedIntervalNames = binIntervals.intervalNames
        File? binnedIntervalsZip = binIntervals.binnedIntervalsZip
    }

    parameter_meta{
        inputRefTarball: {help: "Reference fasta tarball"}
        variantFOF: {help: "Variant file of files. Contains compressed vcf file names"}
        caller: {help: "Variant caller used"}
        intervalBED: {help: "BED file to separate variants into for merging"}
    }
}

task binIntervals {
    input {
        File bedFile
        Int nucleotidesPerBin=85000000
        String oDir="intervals"

        String docker="cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 50,
                    "cpuCount" : 2,
                    "diskGiB" : 2,
                    "runtimeMinutes": 60,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "SSD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    command {
        set -euxo pipefail
        bin_intervals.py \
        --bedfile ~{bedFile} \
        --nucleotides_per_bin ~{nucleotidesPerBin} \
        -odir ~{oDir}

        cat ~{oDir}/binned_intervals_FOF.txt | sed 's!.*/!!' | sed 's/.bed//g' > ~{oDir}/binned_intervals.txt

        zip ~{oDir}/binned_intervals.zip ~{oDir}/*.bed
    }

    output{
        File intervalNames = "~{oDir}/binned_intervals.txt"
        Array[File] binnedIntervals = glob("~{oDir}/*.bed")
        File binnedIntervalsZip = "~{oDir}/binned_intervals.zip"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])} ~{select_first([runtimeAttributesOverride.diskType, defaultRuntimeAttributes.diskType])}"
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

task splitVariantByIntervalBins {
    input {
        Variant variant
        Array[File] intervalBins
        String docker = "cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 25,
                    "cpuCount" : 16,
                    "diskGiB" : 0,
                    "runtimeMinutes": 180,
                    "maxPreemptAttempts": 3,
                    "hpcQueue": "norm",
                    "diskType" : "HDD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(size(variant.variant,  "GB")*2.5) + ceil(size(intervalBins,  "GB")) + 10 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    String autoSampleID = basename(basename(variant.variant, ".g.vcf.gz"), ".genome.vcf.gz")
    # String oDir = "splitVariants"
    
    command <<<
        set -euxo pipefail
        
        cat ~{write_lines(intervalBins)} | \
        parallel -j ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} 'intervalName=$(basename "{.}");  \
            bcftools filter --regions-file {} -o ~{autoSampleID}.${intervalName}.g.vcf.gz ~{variant.variant}; \
            printf "%s\t%s\n" $intervalName ~{autoSampleID}.$intervalName.g.vcf.gz' > variants_interval_files.txt
    >>>

    
    output {
        Map[String, File] variantIntervalMap = read_map("variants_interval_files.txt")
        Array[File] variantInterval = glob("*.g.vcf.gz")
        
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

task glnexus {
    input {
        String caller
        Array[File] variants
        File? intervalBED # Interval bed file with 3 columns: chromosome, start, stop
        String? intervalName
        

        # Runtime Inputs
        String docker="cgrlab/glnexus:v1.4.1"
        String glnexusConfigOption=""
        File? glnexusConfigFile
        String dbDir = "GLNEXUS"
        RuntimeAttributes? runtimeAttributes
    }
    RuntimeAttributes defaultRuntimeAttributes = {
                    "memoryGiB" : 150,
                    "cpuCount" : 32,
                    "diskGiB" : 0,
                    "runtimeMinutes": 600,
                    "maxPreemptAttempts": 0,
                    "hpcQueue": "norm",
                    "diskType" : "HDD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    File variantList = write_lines(variants)
    # String glnexusConfig = if glnexusConfigOption != "" then glnexusConfigOption else basename(glnexusConfigFile)
    String glnexusConfig = if defined(glnexusConfigFile) then basename(select_first([glnexusConfigFile])) else glnexusConfigOption
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(8 * size(variants,  "GB")) + 50 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])

    String oDir = "joint_genotype/~{caller}"
    String oFile = if defined(intervalName) then "~{oDir}/~{caller}.~{intervalName}.glnexus.bcf" else "~{oDir}/~{caller}.glnexus.bcf"

    Int glnexusRAM=floor(select_first([runtimeAttributesOverride.memoryGiB, defaultRuntimeAttributes.memoryGiB])*0.95)
    
    command {
        set -euox pipefail

        mkdir -p ~{oDir}
        ~{"ln -s " + glnexusConfigFile + " ./" }

        glnexus_cli \
        ~{"--bed " + intervalBED} \
        --dir ~{dbDir} \
        --config ~{glnexusConfig} \
        -m ~{glnexusRAM} \
        -t ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} \
        ~{sep=" " variants} \
        > ~{oFile}
        # --list ~{variantList} \

        du --max-depth=1 --human-readable 
        rm -r ~{dbDir}
    }

    output {
        File glnexusBCF="~{oFile}"
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

task indexBCF {
    input {
        File bcf
        String oDir
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
                    "diskType" : "HDD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(size(bcf,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    String basename = basename(bcf)
    
    String bcfIndexFilename = "~{oDir}/~{basename}.csi"
    
    
    command {
        set -euxo pipefail
        mkdir -p ~{oDir} 
        bcftools index --threads ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} -o ~{bcfIndexFilename} ~{bcf} 
        cp ~{bcf} ~{oDir}/~{basename}
    }
    output {
        File bcfOriginal = "~{oDir}/~{basename}"
        File bcfIndex = "~{bcfIndexFilename}"
        # File bcfZip = "~{basename}.zip"
        
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

task concatVariantIntervals {
    input {
        BcfAndIndexArray bcfsAndIndices
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
                    "diskType" : "HDD"
    }
    RuntimeAttributes runtimeAttributesOverride = select_first([runtimeAttributes, defaultRuntimeAttributes])
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(2.0 * size(bcfsAndIndices.bcf,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) 
    String oDir="joint_genotype/~{caller}"
    String oFile = "~{oDir}/~{caller}.concat.bcf"

    command {
        set -euxo pipefail
        mkdir -p ~{oDir} 
        bcftools concat  --threads ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} -Ou -a ~{sep=' ' bcfsAndIndices.bcf} > ~{oFile}
    }
    output {
        File concatVariants = "~{oFile}"
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

task compressBCF {
    input {
        File bcf
        # String caller
        String docker = "cgrlab/gencompass@sha256:ab81377e6c793d8bbf5b44bc1ea4e5a4b1bec943b7fb0691fd02ca9cfc64c21a"
        String oDir="joint_genotype"
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
    Int autoDiskGB = if select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB]) < 1 then ceil(2.0 * size(bcf,  "GB")) + 65 else select_first([runtimeAttributesOverride.diskGiB, defaultRuntimeAttributes.diskGiB])
    String basename = basename(bcf)
    # String oDir="~{caller}/genotyped"

    
    command {
        set -euxo pipefail
        mkdir -p ~{oDir}
        cd ~{oDir} 
        bcftools view --min-ac 1 --threads ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} ~{bcf} | \
        bgzip -c --threads ~{select_first([runtimeAttributesOverride.cpuCount, defaultRuntimeAttributes.cpuCount])} > ~{basename}.gz; \
        tabix -p vcf ~{basename}.gz
    }
    output {
        File bcfCompressed = "~{oDir}/~{basename}.gz"
        File bcfIndex = "~{oDir}/~{basename}.gz.tbi"
        # File bcfZip = "~{basename}.zip"
        
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

struct BcfAndIndexArray {
    Array[File] bcf
    Array[File] bcfIndex
}

struct BCFandIndex {
    File bcf
    File bcfIndex
}

struct Variant {
    File variant
    File variantIndex
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


