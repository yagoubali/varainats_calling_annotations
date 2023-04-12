#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.fastqDir="/media/yagoubali/bioinfo3/Project_Paul/SRR22797223"
params.fastqDir="/media/yagoubali/bioinfo3/db/ibitayo/fastq"
params.outDir="testing"
raw_reads = params.fastqDir
out_dir = file(params.outDir,  mode: "copy")

out_dir.mkdir()


process runFastQC{
    cpus { 2 }
    publishDir "${out_dir}/qc/raw/${pair_id}", mode: 'copy', overwrite:false

    input:
    //Channel.fromFilePairs("${raw_reads}/*R[1,2].001.fastq.gz", type: 'file')
    tuple val(pair_id), path(reads) 
   

    output:
        //file("${sample}_fastqc/*.zip") //into fastqc_files
        //file("test.txt")
        file("${pair_id}_fastqc/*.zip")
        //path("${publishDir}/**/*R[1,2]_001.fastq.gz.zip")

    
    script:
    """
     echo "${pair_id} ... ${reads[0]} ${reads[1]}" > test.txt
     mkdir ${pair_id}_fastqc
     echo "${pair_id}" > ${pair_id}_fastqc/${reads[0]}.zip
     echo "${pair_id}" > ${pair_id}_fastqc/${reads[1]}.zip
    """
    /*
     fastqc --outdir ${pair_id}_fastqc \
      ${params.fastqDir}/${reads[0]} \
      ${params.fastqDir}/${reads[1]
     */ 
}

workflow {
    read_pair = Channel.fromFilePairs("${raw_reads}/*R[1,2]_001.fastq.gz", type: 'file')
    println read_pair
    pair_id_qc_channel = runFastQC(read_pair).map{T->[T[0],T[1]]}
    println pair_id_qc_channel.view()


   
}