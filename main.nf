#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.fastqDir="/media/yagoubali/bioinfo3/Project_Paul/SRR22797223"
params.fastqDir="/media/yagoubali/bioinfo3/Project_Paul/test_data"
params.reference_baseFolder="/media/yagoubali/bioinfo3/Project_Paul/broadinstitute/bundle"
referenceFile="${params.reference_baseFolder}/Reference/Homo_sapiens_assembly38.fasta"
params.gatkPath="/media/yagoubali/bioinfo3/Project_Paul/gatk-4.4.0.0/gatk"
params.outDir="testing"
params.pl= "illumina"
raw_reads = params.fastqDir
out_dir = file(params.outDir,  mode: "copy")

reference=channel.value(referenceFile)
out_dir.mkdir()


process runFastQC{
    cpus { 2 }
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
    publishDir "${out_dir}/markDuplicates", mode: 'copy', overwrite:false
    
    input:
        tuple val(pair_id), path(reads) 

    output:
        path ("*dedup_metrics.txt")
        path("*sorted_dedup.bam")


    script:
    ///home/yagoubali/anaconda3/bin/gatk  --java-options "-Djava.io.tmpdir=tmpdir" 
    //MarkDuplicatesSpark       -I SRR14509522.sorted.bam    -M dedup_metrics.txt    -O sorted_dedup.bam  --remove-sequencing-duplicates

    """
    mkdir -p ${pair_id}
    ${params.gatkPath}  \
     --java-options "-Djava.io.tmpdir=${pair_id}" \
	 MarkDuplicatesSpark \
	-I ${out_dir}/mapping/${pair_id}/${reads[0]} \
	-M ${pair_id}_dedup_metrics.txt \
	-O ${pair_id}_sorted_dedup.bam 
    """ 
    // rm -r ${pair_id}
}

process run_BaseRecalibrator{
     cpus { 2 }
    memory '2 GB'
    //publishDir "${out_dir}/markDuplicates", mode: 'copy', overwrite:false
    
    input:
        tuple val(pair_id), path(reads) 

}
workflow {
    read_pair = Channel.fromFilePairs("${raw_reads}/**/*R[1,2].fastq.gz", type: 'file')
    //println read_pair
    // --> pair_id_qc_channel = runFastQC(read_pair).map{T->[T[0],T[1]]}
    // ---> println pair_id_qc_channel.view()
    qcFiles=runFastQC(read_pair)
    multiqc(qcFiles)
    //println read_pair.view()
    read_pair2 = Channel.fromFilePairs("${raw_reads}/**/*R[1,2].fastq.gz", type: 'file')
    genomeSize_heterozygosity(read_pair)
    //println read_pair.view()
    //pair_id_mapped_channel = 
    sorted_bam_ch=run_bwa(read_pair2).map{T->[T[0], T[1]]}
    //.fromPath( "${out_dir}/mapping/**/*sorted.bam")
    // println pair_id_mapped_channel.view()
     println sorted_bam_ch.view()
     run_markDuplicatesSpark(sorted_bam_ch)

///home/yagoubali/anaconda3/bin/gatk 
}