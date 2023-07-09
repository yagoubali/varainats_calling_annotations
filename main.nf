#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// changes lines 5-10
params.fastqDir="/media/yagoubali/bioinfo3/Project_Paul/test_data"
params.reference_baseFolder="/media/yagoubali/bioinfo3/Project_Paul/broadinstitute/bundle"
referenceFile="${params.reference_baseFolder}/Reference/Homo_sapiens_assembly38.fasta"
params.outDir="testing"
params.pl= "illumina"
params.maxForks=2

// No need to make any changes in the below lines
raw_reads = params.fastqDir
out_dir = file(params.outDir,  mode: "copy")
reference=channel.value(referenceFile)
out_dir.mkdir()
params.gatkPath="~/bin/gatk-4.4.0.0/gatk"
params.snpEFF="~/bin/snpEff/snpEff.jar"
process runFastQC{
    cpus { 2 }
    maxForks 2
    container "yagoubali/fastqc:0.11.9"
    publishDir "${out_dir}/qc/raw/${pair_id}", mode: 'copy', overwrite:false

    input:
    //Channel.fromFilePairs("${raw_reads}/*R[1,2].001.fastq.gz", type: 'file')
    tuple val(pair_id), path(reads) 
   

    output:
        path("${pair_id}_fastqc/*.zip") //into fastqc_files
    
    script:

    """
     mkdir ${pair_id}_fastqc
     fastqc --outdir ${pair_id}_fastqc  \
     ${params.fastqDir}/${pair_id}/${reads[0]} \
     ${params.fastqDir}/${pair_id}/${reads[1]}

    """

}

process multiqc{
    container 'yagoubali/multiqc'
    publishDir "${out_dir}/qc/raw", mode: 'copy', overwrite: false

    input:
        path('*.zip')

    output:
        file('multiqc_report.html')

    """
    multiqc ${out_dir}/qc/raw
    """
}

process genomeSize_heterozygosity {
    cpus { 2 }
    memory '2 GB'
    maxForks 2
    container 'yagoubali/jellyfish:2.2.10'

    publishDir "${out_dir}/genomeSize_heterozygosity/${pair_id}", mode: 'copy', overwrite:false

    input:
    //Channel.fromFilePairs("${raw_reads}/*R[1,2].001.fastq.gz", type: 'file')
    tuple val(pair_id), path(reads) 
     output:
     path("log.log")
     path("${pair_id}_pe12")
     path("${pair_id}_mer_counts_dumps.fa")
     path("${pair_id}_hito")
     path("figures", type: 'dir')
     //path("${pair_id}.fastq")
    
    script:
    """
    echo "jellyfish count -t 30 -C -m 21 -s 5G -o ${pair_id}_pe12 <(zcat ${params.fastqDir}/${pair_id}/${reads[0]}) <(zcat ${params.fastqDir}/${pair_id}/${reads[1]})" > log.log
    jellyfish count -t 2 -C -m 21 -s 100M -o ${pair_id}_pe12  <(zcat ${params.fastqDir}/${pair_id}/${reads[0]}) <(zcat ${params.fastqDir}/${pair_id}/${reads[1]})
    jellyfish dump ${pair_id}_pe12 > ${pair_id}_mer_counts_dumps.fa
    jellyfish histo -o ${pair_id}_hito ${pair_id}_pe12
    genomescope.R -i ${pair_id}_hito  -o figures -k 21

    """

}

process run_bwa {
    cpus { 2 }
    memory '2 GB'
    maxForks 2
    container 'yagoubali/bwa'

    publishDir "${out_dir}/mapping/${pair_id}", mode: 'copy', overwrite:false

    input:
    //Channel.fromFilePairs("${raw_reads}/*R[1,2].001.fastq.gz", type: 'file')
    tuple val(pair_id), path(reads) 
    output:
    tuple val(pair_id), path ("*.sorted.bam"), path ("*.bam"), path ("*.flagstat")

    script:
    readGroup="@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:${params.pl}\\tSM:${pair_id}"
    """
    bwa mem -t 2 -aM -R \"${readGroup}\" ${referenceFile}   \
    ${params.fastqDir}/${pair_id}/${reads[0]}   \
    ${params.fastqDir}/${pair_id}/${reads[1]} > ${pair_id}.unsorted.raw.sam
    samtools view -S -b ${pair_id}.unsorted.raw.sam > ${pair_id}.bam
    samtools sort -o ${pair_id}.sorted.bam ${pair_id}.bam 
    samtools flagstat ${pair_id}.bam > ${pair_id}.flagstat

    """
}

process run_markDuplicatesSpark {
    cpus { 2 }
    memory '2 GB'
    maxForks 2
    container 'yagoubali/gatk:4.4.0.0'

    publishDir "${out_dir}/markDuplicates", mode: 'copy', overwrite:false
    
    input:
        tuple val(pair_id), path(reads) 

    output:
         tuple val(pair_id), path("*dedup_metrics.txt"), path("*sorted_dedup.bam")


    script:
    """
    mkdir -p ${pair_id}
    ${params.gatkPath}  \
     --java-options "-Djava.io.tmpdir=${pair_id}" \
	 MarkDuplicatesSpark \
	-I ${out_dir}/mapping/${pair_id}/${reads[0]} \
	-M ${pair_id}_dedup_metrics.txt \
	-O ${pair_id}_sorted_dedup.bam 
    """ 
}

process run_BaseRecalibrator_ApplyBQSR{
    cpus { 2 }
    memory '2 GB'
    maxForks 2
    container 'yagoubali/gatk:4.4.0.0'

    publishDir "${out_dir}/BaseRecalibrator_ApplyBQSR", mode: 'copy', overwrite:false
    
    input:
        tuple val(pair_id), path(reads) 
    output:
        tuple val(pair_id), path("*_recal_data.table"), path("*_dup.bqsr.bai"), path("*_dup.bqsr.bam"), path("*vcf.gz"), path("*.g.vcf.gz.tbi") 
    
    script:    

    """
      mkdir -p ${pair_id}
    ${params.gatkPath}  --java-options "-Djava.io.tmpdir=${pair_id}" \
     BaseRecalibrator \
    -I ${out_dir}/markDuplicates/${pair_id}_sorted_dedup.bam \
    -R ${params.reference_baseFolder}/Reference/Homo_sapiens_assembly38.fasta \
    --known-sites ${params.reference_baseFolder}/hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
    -O ${pair_id}_recal_data.table


    ${params.gatkPath}    --java-options "-Djava.io.tmpdir=${pair_id}" \
    ApplyBQSR \
    -I ${out_dir}/markDuplicates/${pair_id}_sorted_dedup.bam \
    -R ${params.reference_baseFolder}/Reference/Homo_sapiens_assembly38.fasta \
    --bqsr-recal-file ${pair_id}_recal_data.table \
    -O ${pair_id}_dup.bqsr.bam

    ${params.gatkPath}  --java-options "-Djava.io.tmpdir=${pair_id}" \
     HaplotypeCaller \
     -I ${pair_id}_dup.bqsr.bam \
     -R ${params.reference_baseFolder}/Reference/Homo_sapiens_assembly38.fasta \
     -ERC GVCF \
    -O ${pair_id}.g.vcf.gz
    """  

}

process run_prepareVCF{
    cpus { 2 }
    memory '2 GB'
    maxForks 2
    container 'yagoubali/gatk:4.4.0.0'

    publishDir "${out_dir}/final", mode: 'copy', overwrite:false

    input:
    val x
        //tuple val(pair_id), path(reads)
       
    output:
        tuple path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"), path("output.vcf.gz"), path("output.vcf.gz.tbi"), path("snps.recal"), path ("snps.tranches"), path("output.vqsr.vcf"), path("output.vqsr.vcf.idx")

    
    script:    

    """
    mkdir ToCombine
     ls ${out_dir}/BaseRecalibrator_ApplyBQSR/*g.vcf.gz > ToCombine.list
    ${params.gatkPath}  --java-options "-Djava.io.tmpdir=ToCombine" \
     CombineGVCFs \
    -R ${params.reference_baseFolder}/Reference/Homo_sapiens_assembly38.fasta \
    -V ToCombine.list \
    -O combined.g.vcf.gz


    ${params.gatkPath}  --java-options "-Djava.io.tmpdir=ToCombine" \
     GenotypeGVCFs \
    -R ${params.reference_baseFolder}/Reference/Homo_sapiens_assembly38.fasta \
    -V combined.g.vcf.gz \
    -O output.vcf.gz

    ${params.gatkPath}  --java-options "-Djava.io.tmpdir=ToCombine" \
     VariantRecalibrator  -V output.vcf.gz  \
    --trust-all-polymorphic \
    -mode SNP \
    --max-gaussians 6 \
    --resource:dbsnp,known=true,training=true,truth=true,prior=7  ${params.reference_baseFolder}/hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -O snps.recal \
    --tranches-file snps.tranches

    ${params.gatkPath}  --java-options "-Djava.io.tmpdir=ToCombine" \
    ApplyVQSR \
    -R ${params.reference_baseFolder}/Reference/Homo_sapiens_assembly38.fasta \
    -V output.vcf.gz \
    -O output.vqsr.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file snps.tranches \
    --recal-file snps.recal \
    -mode SNP
    """  

}


process run_snpEFF{
    cpus { 2 }
    memory '2 GB'
    maxForks 2
    container 'yagoubali/snpeff'


    publishDir "${out_dir}/final", mode: 'copy', overwrite:false

    input:
    val x
        //tuple val(pair_id), path(reads)
       
    output:
        tuple path("output_annotated.vcf"), path("snpEff_summary.html"),path("snpEff_genes.txt")


    
    script:    

    """
    java -Xmx8g -jar  ${params.snpEFF} GRCh38.105 \
    ${out_dir}/final/output.vqsr.vcf>output_annotated.vcf
    """
}

workflow {
    read_pair = Channel.fromFilePairs("${raw_reads}/**/*R[1,2].fastq.gz", type: 'file')
    qcFiles=runFastQC(read_pair)
    multiqc(qcFiles)
    read_pair2 = Channel.fromFilePairs("${raw_reads}/**/*R[1,2].fastq.gz", type: 'file')
    genomeSize_heterozygosity(read_pair)
    sorted_bam_ch=run_bwa(read_pair2).map{T->[T[0], T[1]]}
    BaseRecalibrator_ApplyBQSR_ch=run_markDuplicatesSpark(sorted_bam_ch).map{T->[T[0],T[2]]}
    prepare_vcf_ch=run_BaseRecalibrator_ApplyBQSR(BaseRecalibrator_ApplyBQSR_ch)
    annotation_ch=run_prepareVCF(prepare_vcf_ch.collect())
    snpEFF_ch=run_snpEFF(annotation_ch.collect())
    //println snpEFF_ch.view()
    
}
