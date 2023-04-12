#!/usr/bin/env bash

##prefetch SRP319461
fastqDir="/media/yagoubali/bioinfo3/Project_Paul/data";
subDirs=($(ls ${fastqDir}))
for  Dir in "${subDirs[@]}"; do
cd ${fastqDir}/${Dir};
ls ${Dir}.sra
echo " Starting fastq-dump .....";
fastq-dump --gzip --split-files --readids  ${Dir}.sra
echo " Starting renaming output file to R1, R2 .....";
rename 's/_1/_R1/' *.fastq.gz
rename 's/_2/_R2/' *.fastq.gz
cd ..;
#
#ls ${fastqDir}/${Dir}
 done 