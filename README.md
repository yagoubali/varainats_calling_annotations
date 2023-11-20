**README for NGS Data Analysis Pipeline
Background**
This pipeline is designed to provide a simplified and easy-to-use tool for genomic analysis. It covers different steps of the workflow for the analysis of NGS data including quality control and genome size checks, variant calling and annotation. It is built using Nextflow which is a workflow language designed at the Centre for Genomic Regulation, Barcelona. The knowledge of Nextflow is not necessary to use this pipeline, as the heavy-lifting job of putting together the script has been done by the designers of this pipeline. Our pipeline also allows for executions to resume at the point where the workflow was stopped in the case of additional inputs or configuration changes rather than starting afresh. 

Docker has also been used as a containerizing tool to assemble all the software and dependencies needed for the analysis, thereby, ensuring an easy installation and execution of the pipeline. As such, the user only needs to make an installation of this pipeline and then good to go on with the analysis because it automatically installs all the software dependencies.  

In addition, time-consuming steps, such as indexing of a reference sequence, are only performed once when executing the workflow multiple times using the same working directory. 

Nextflow and the Docker software can be run on Linux, Mac OS, and Windows supported through WSL2.

Java 17 was required for the gatk4 to work effectively. 

**Quality Control Checks**
The pipeline is designed to be able to execute and generate quality reports using both fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and Multiqc (https://doi.org/10.1093/bioinformatics/btw354) workflows effectively.   

The output is directly piped into samtools, converted into BAM format, sorted, and then indexed.
Genome size heterozygosity
We implemented scripts for determining the genome characteristics such as heterozygosity, repeat content, and size from sequencing reads using Jellyfish and GenomeScope. 

**Variant Calling**
The variant calling step is carried out using the Genome Analysis Toolkit (GAKT) best practices workflow. This is done together with samtools, bwa-mem. 

**Variant Annotation**
The annotation of the variants is done using SnpEff. 
