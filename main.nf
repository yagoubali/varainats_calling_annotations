#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.fastqDir="/media/yagoubali/bioinfo3/Project_Paul/SRR22797223"
params.fastqDir="/media/yagoubali/bioinfo3/Project_Paul/test_data"
params.reference_baseFolder="/media/yagoubali/bioinfo3/Project_Paul/broadinstitute/bundle"
referenceFile="${params.reference_baseFolder}/Reference/Homo_sapiens_assembly38.fasta"
params.gatkPath="/media/yagoubali/bioinfo3/Project_Paul/gatk-4.4.0.0/gatk"
params.outDir="testing"
params.pl= "illumina"
params.maxForks=2
raw_reads = params.fastqDir
out_dir = file(params.outDir,  mode: "copy")

reference=channel.value(referenceFile)
out_dir.mkdir()


process runFastQC{
    cpus { 2 }
    maxForks 2
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

    publishDir "${out_dir}/BaseRecalibrator_ApplyBQSR", mode: 'copy', overwrite:false
    
    input:
        tuple val(pair_id), path(reads) 
    output:
        tuple val(pair_id), path("*_recal_data.table"), path("*_dup.bqsr.bai"), path("*_dup.bqsr.bam"), path("*vcf.gz") 
    
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

process run_CombineGVCFs_GenotypeGVCFs{
    cpus { 2 }
    memory '2 GB'
    maxForks 2

    publishDir "${out_dir}/HaplotypeCaller", mode: 'copy', overwrite:false
    
    input:
        tuple val(pair_id), path(reads) 
    output:
        tuple val(pair_id), path("*vcf.gz")  
    
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
    """  

}
workflow {
    read_pair = Channel.fromFilePairs("${raw_reads}/**/*R[1,2].fastq.gz", type: 'file')
    //println read_pair
    // --> pair_id_qc_channel = runFastQC(read_pair).map{T->[T[0],T[1]]}
    // ---> println pair_id_qc_channel.view()
    qcFiles=runFastQC(read_pair)
    multiqc(qcFiles)
    read_pair2 = Channel.fromFilePairs("${raw_reads}/**/*R[1,2].fastq.gz", type: 'file')
    genomeSize_heterozygosity(read_pair)
    sorted_bam_ch=run_bwa(read_pair2).map{T->[T[0], T[1]]}
    BaseRecalibrator_ApplyBQSR_ch=run_markDuplicatesSpark(sorted_bam_ch).map{T->[T[0],T[2]]}
    run_BaseRecalibrator_ApplyBQSR(BaseRecalibrator_ApplyBQSR_ch)

    println BaseRecalibrator_ApplyBQSR_ch.view()

}