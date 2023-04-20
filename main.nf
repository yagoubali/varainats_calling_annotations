#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.fastqDir="/media/yagoubali/bioinfo3/Project_Paul/SRR22797223"
params.fastqDir="/media/yagoubali/bioinfo3/Project_Paul/test_data"
params.reference="/media/yagoubali/bioinfo3/Project_Paul/Reference_data/resource-GRCh38/hs38DH-extra.fa"
params.outDir="testing"
raw_reads = params.fastqDir
out_dir = file(params.outDir,  mode: "copy")

reference=channel.value(params.reference)
out_dir.mkdir()


process runFastQC{
    cpus { 2 }
    publishDir "${out_dir}/qc/raw/${pair_id}", mode: 'copy', overwrite:false

    input:
    //Channel.fromFilePairs("${raw_reads}/*R[1,2].001.fastq.gz", type: 'file')
    tuple val(pair_id), path(reads) 
   

    output:
        path("${pair_id}_fastqc/*.zip") //into fastqc_files
        //file("test.txt")
        // tuple val(pair_id), file("${pair_id}_fastqc/*.zip")
        //path("${publishDir}/**/*R[1,2]_001.fastq.gz.zip")

    
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
    path("*.bam")
    path("*.flagstat")

    script:
    """
    bwa mem -t 2 ${params.reference}   \
    ${params.fastqDir}/${pair_id}/${reads[0]}   \
    ${params.fastqDir}/${pair_id}/${reads[1]} > ${pair_id}.unsorted.raw.sam
    samtools view -S -b ${pair_id}.unsorted.raw.sam > ${pair_id}.bam
    samtools sort -o ${pair_id}.sorted.bam ${pair_id}.bam 
    samtools flagstat ${pair_id}.bam > ${pair_id}.flagstat

    """
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
    println read_pair.view()
    run_bwa(read_pair2)

}