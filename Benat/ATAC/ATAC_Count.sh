#!/bin/bash
######## Variable names ##############
### Files array
# list fastq files | look for "fastq.gz" format | remove the last x (16) characters
FILES=($(ls /home/bariceta/RINEY/ATAC/BAM/ | grep ".sort.rmdup.rmblackls.rmchr.bam$" | sed -r 's/.{31}$//' | uniq ))
### root folder to work
BaseFolder=/home/bariceta/RINEY/ATAC/
### Is it paired end sequencing? 0 = No 1 = Yes
PairedEnd=1

### count matrix path and name
# store as .Rdata file
countMat="/home/bariceta/RINEY/ATAC/Rsubread_Counts_ATAC_Study1.Rdata"
### Reference genome folder
Genome="/datos/intercambio/GRCh38_v103/"
### path and folder where scripts are placed
ScriptsFolder="/home/bariceta/scripts_new/ATAC/"
### export file names to .txt 
echo "${FILES[@]}" > ${BaseFolder}files_count.txt
##########################




# Count
sbatch --output=${BaseFolder}Count.log --export=ALL,FILES="${BaseFolder}files_count.txt",BaseFolder=$BaseFolder,PairedEnd=$PairedEnd,countMat=$countMat,Genome=$Genome  ${ScriptsFolder}SlurmBatch_R.sbs
