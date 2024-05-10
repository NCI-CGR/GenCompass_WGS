version 1.0

workflow Annotate {
    input{
        String name
        File vcf
        Array[File] annovarDatabases
        Array[String] annovarProtocols = ["refGene","gnomad312_genome","clinvar_20230416","cytoBand"]
        Array[String] annovarOperations = ["g","f","f","r"]

        Array[File] annovarForIntervarDatabases
        Array[String] annovarForIntervarProtocols = ["refGene","esp6500siv2_all","avsnp147","dbnsfp42a","clinvar_20230416","gnomad_genome","dbscsnv11","rmsk","ensGene","knownGene"]
        Array[String] annovarForIntervarOperations = ["g","f","f","f","f","f","f","r","g","g"]

        String snpEffDatabase="GRCh38.105"
        File dbSnpIdTable

        String annovarDocker
        String snpEffDocker
        String gencompassDocker
    }
    VCF dbSnpIDVCF = {
        "vcf": dbSnpIdTable,
        "tbi": dbSnpIdTable+".tbi"
    }
    
    call convertToAnnovar{
    input:
        vcf=vcf,
        docker=annovarDocker
    }

    call annovarAnnotate {
        input:
            avinput = convertToAnnovar.avinput,
            databases = annovarDatabases,
            protocols=annovarProtocols,
            operations=annovarOperations,
            docker=annovarDocker,
            name=name
    }

    call annovarAnnotate as annovarAnnotateForIntervar {
        input:
            avinput = convertToAnnovar.avinput,
            databases = annovarForIntervarDatabases,
            protocols=annovarForIntervarProtocols,
            operations=annovarForIntervarOperations,
            name=name+"_intervar",
            docker=annovarDocker
    }
    call intervarAnnotate {
        input:
            avinput = annovarAnnotateForIntervar.annovarAnnotations,
            # databases = intervarAnnovarDatabases,
            docker = annovarDocker,
            name=name

    }

    call snpEffAnnotate {
        input:
            vcf = vcf,
            name=name,
            database=snpEffDatabase,
            docker=snpEffDocker
    }
    call snpEffAddDbSnpIdField{
        input:
            vcf = snpEffAnnotate.snpEffAnnot,
            docker = snpEffDocker,
            dbSnpIDVCF=dbSnpIDVCF
    }
    call bcftoolsExtractSnpEffAnn {
        input:
            vcf = snpEffAddDbSnpIdField.vcfWithSnpId,
            name=name,
            docker=gencompassDocker
    }
    call parseIntervar {
        input:
            name=name,
            intervar = intervarAnnotate.intervar,
            docker = gencompassDocker
    }

    call parseSnpEff {
        input:
            name=name,
            snpEffTSV=bcftoolsExtractSnpEffAnn.extractedSnpEff,
            docker = gencompassDocker
    }

    call mergeAnnotations{
        input:
            name = name,
            annovar = annovarAnnotate.annovarAnnotations,
            parsedSnpEff=parseSnpEff.parsedSnpEff,
            parsedIntervar=parseIntervar.parsedIntervar,
            docker = gencompassDocker
    }
    output {
        File annovarAnnotations = annovarAnnotate.annovarAnnotations
        File annovarIntervarAnnotations = annovarAnnotateForIntervar.annovarAnnotations
        File intervarAnnotations = intervarAnnotate.intervar
        File snpEffVCF = snpEffAddDbSnpIdField.vcfWithSnpId
        
        File geneAnnotation = mergeAnnotations.geneAnnotation
        File intergenicAnnotation = mergeAnnotations.intergenicAnnotation 
        File intronicAnnotation = mergeAnnotations.intronicAnnotation 
        File ncRNAAnnotation = mergeAnnotations.ncRNAAnnotation
        File fullMergedAnnotation = mergeAnnotations.fullMergedAnnotation

    }
}

task convertToAnnovar{
    input{
        File vcf
        
        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        Int gbRAM = 10

    }
    Int diskSizeGiB = ceil(size(vcf, "GiB") * 2.5)
    String avinputFilename = basename(vcf, ".vcf") + ".avinput"
    command <<<
        convert2annovar.pl -format vcf4 --keepindelref ~{vcf} > ~{avinputFilename}
    >>>
    output{
        File avinput = "~{avinputFilename}"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{diskSizeGiB} SSD"
        cpu : 2
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts

    }
}
task annovarAnnotate {
    input{
        File avinput
        Array[File] databases
        String buildVersion = "hg38"
        String name
        Array[String] protocols = ["refGene","gnomad312_genome","clinvar_20221231","cytoBand"]
        Array[String] operations = ["g","f","f","r"]
        String naString = "."

        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        Int gbRAM = 10

    }
    Int diskSizeGiB = ceil(size(avinput, "GiB") * 2.5)
    command <<<
        set -euox pipefail
        mkdir -p annovar_db
        ln -s ~{sep =" " databases} ./annovar_db/ 

        table_annovar.pl \
        --buildver ~{buildVersion} \
        --out ~{name}_annovar \
        --remove \
        --verbose \
        --protocol ~{sep=',' protocols} \
        --operation ~{sep=',' operations} \
        --nastring ~{naString} \
        ~{avinput} \
        ./annovar_db 

        sleep 5
    >>>
    output{
        File annovarAnnotations = "~{name}_annovar.~{buildVersion}_multianno.txt"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{diskSizeGiB} SSD"
        cpu : 2
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts

    }
}

task intervarAnnotate {
    input{
        File avinput
        # Array[File] databases
        String buildVersion = "hg38"
        String name


        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        Int gbRAM = 10

    }
    Int diskSizeGiB = ceil(size(avinput, "GiB") * 2.5)
    command <<<
        set -euo pipefail

        ln -s ~{avinput} ~{name}.~{buildVersion}_multianno.txt

        Intervar.py \
            -i ~{avinput} \
            --skip_annovar \
            -b ~{buildVersion} \
            -o ~{name}
        
        sleep 5
    >>>
    output{
        File intervar = "~{name}.hg38_multianno.txt.intervar"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{diskSizeGiB} SSD"
        cpu : 2
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts

    }
}

task snpEffAnnotate {
    input{
        File vcf
        File? config
        String database="GRCh38.105"
        String name
        Array[String] annotationParameters = [
            "-no-intergenic",
            "-no INTRAGENIC",
            "-no-downstream",
            "-no-upstream", 
            "-no-utr", 
            "-no PROTEIN_STRUCTURAL_INTERACTION_LOCUS", 
            "-no NEXT_PROT", 
            "-lof", 
            "-noStats" ]

        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        Int gbRAM = 100
    }
    Int diskSizeGiB = ceil(size(vcf, "GiB")*2.5)
    command <<<
    set -euox pipefail
    snpEff \
    ~{"-config " + config} \
    -verbose \
    ~{sep = " " annotationParameters} \
    ~{database} \
    ~{vcf} > ~{name}.snpEff.vcf
    >>>

    output {
        File snpEffAnnot = "~{name}.snpEff.vcf"

    }
    runtime {
        docker : docker
        disks : "local-disk ~{diskSizeGiB} SSD"
        cpu : 2
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts

    }
}

task snpEffAddDbSnpIdField{
    input {
        File vcf
        VCF dbSnpIDVCF

        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        Int gbRAM = 10
    }
    String vcfBasename = basename(vcf, ".vcf")
    Int diskSizeGiB = ceil(size(vcf, "GiB") * 2 + size(dbSnpIDVCF.vcf, "GiB") * 1.2)
    command <<<
    set -euox pipefail
    SnpSift annotate -id ~{dbSnpIDVCF.vcf} ~{vcf} > ~{vcfBasename}.dbSNPID.vcf
    >>>
    output {
        File vcfWithSnpId = "~{vcfBasename}.dbSNPID.vcf"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{diskSizeGiB} SSD"
        cpu : 2
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts

    }

}

task bcftoolsExtractSnpEffAnn {
    input {
        File vcf
        String name

        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        Int gbRAM = 10
    }
    Int diskSizeGiB = ceil(size(vcf, "GiB") * 2.5 )
    command <<<
    set -euox pipefail
    bcftools query  --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/set\t%INFO/ANN\n' \
        ~{vcf} > ~{name}.snpEff.tsv
    >>>
    output {
        File extractedSnpEff = "~{name}.snpEff.tsv"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{diskSizeGiB} SSD"
        cpu : 2
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts

    }

}

task parseIntervar {
    input {
        File intervar
        String name
        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        Int gbRAM = 10
    }
    Int diskSizeGiB = ceil(size(intervar, "GiB") * 2.5 )
    command <<<
    set -euox pipefail
    parse_intervar.py \
    --intervar ~{intervar} \
    --output ~{name}.intervar_parsed.tsv
    >>>
    output {
        File parsedIntervar = "~{name}.intervar_parsed.tsv"

    }
    runtime {
        docker : docker
        disks : "local-disk ~{diskSizeGiB} SSD"
        cpu : 2
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts

    }

}
task parseSnpEff {
    input {
        File snpEffTSV
        String name
        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        Int gbRAM = 10
    }
    Int diskSizeGiB = ceil(size(snpEffTSV, "GiB") * 2.5 )
    command <<<
        set -euox pipefail
        parse_snpeff.py \
        --include-header \
        --write-unannotated \
        ~{snpEffTSV}  \
        ~{name}.snpEff_parsed.tsv
    >>>
    output {
        File parsedSnpEff = "~{name}.snpEff_parsed.tsv"
    }
    runtime {
        docker : docker
        disks : "local-disk ~{diskSizeGiB} SSD"
        cpu : 2
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts

    }

}

task mergeAnnotations {
    input {
        File annovar
        File parsedIntervar
        File parsedSnpEff
        String name
        String oDir = "annotations"
        String docker
        Int runtimeMinutes = 60
        Int maxPreemptAttempts = 3
        Int gbRAM = 10
    }
    Int diskSizeGiB = ceil(size(annovar, "GiB") * 2.5 + size(parsedIntervar, "GiB") * 2.5 + size(parsedSnpEff, "GiB") * 2.5)
    command <<<
        set -euox pipefail
        mkdir -p ~{oDir}
        merge_annotations.py \
        --annovar ~{annovar} \
        --intervar ~{parsedIntervar} \
        --snpeff ~{parsedSnpEff} \
        --name ~{name} \
        -odir ~{oDir}
    >>>
    output {
        File geneAnnotation = "~{oDir}/~{name}_annotated.gene.tsv"
        File intergenicAnnotation = "~{oDir}/~{name}_annotated.intergenic.tsv" 
        File intronicAnnotation = "~{oDir}/~{name}_annotated.intronic.tsv" 
        File ncRNAAnnotation = "~{oDir}/~{name}_annotated.ncRNA.tsv" 
        File fullMergedAnnotation = "~{oDir}/~{name}_merged_annotations.tsv" 
    }
    runtime {
        docker : docker
        disks : "local-disk ~{diskSizeGiB} SSD"
        cpu : 2
        memory : "~{gbRAM} GiB"
        hpcMemory : gbRAM
        hpcQueue : "norm"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts

    }
}

struct VCF {
    File vcf
    File tbi
}