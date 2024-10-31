########################################################################
################# ATACseq preprocessing script ##########################
########################################################################

#!/bin/bash
######## Variable names ##############
### Files array
# list fastq files | look for "fastq.gz" format | remove the last x (16) characters
FILES=($(ls /home/bariceta/RINEY/ATAC/fastq/ | grep "fastq.gz" | sed -r 's/.{16}$//' | uniq ))
### root folder to work
BaseFolder=/home/bariceta/RINEY/ATAC/
### Is it paired end sequencing? 0 = No 1 = Yes
PairedEnd=1


### Reference genome folder
Genome="/home/bariceta/Alignments/Bowtie2/hg38_Bowtie2/"
### HUMAN
#GenomeIndex=/home/bariceta/Alignments/Bowtie2/hg38_Bowtie2/hg38
#chrOrder=/home/bariceta/Alignments/Bowtie2/hg38_Bowtie2/names.txt
#chrSize=/home/bariceta/Alignments/Bowtie2/hg38_Bowtie2/hg38.chrom.sizes

### path and folder where scripts are placed
ScriptsFolder="/home/bariceta/scripts_new/ATAC/"
### export file names to .txt 
echo "${FILES[@]}" > ${BaseFolder}files.txt
##########################

######## Create folders #########
# fastQC
if [ ! -e ${BaseFolder}fastqc/ ]; then
	mkdir ${BaseFolder}fastqc/
fi
# BAM
if [ ! -e ${BaseFolder}BAM/ ]; then
	mkdir ${BaseFolder}BAM/
fi
# Library metrics
if [ ! -e ${BaseFolder}libMetrics/ ]; then
	mkdir ${BaseFolder}libMetrics/
fi
# summary files
if [ ! -e ${BaseFolder}QC_Align/ ]; then
	mkdir ${BaseFolder}QC_Align/
fi
# Peak Calling
if [ ! -e ${BaseFolder}PeakCalling/ ]; then
	mkdir ${BaseFolder}PeakCalling/
fi
# MultiQC
if [ ! -e ${BaseFolder}multiqc/ ]; then
	mkdir ${BaseFolder}multiqc/
fi


######## Send scripts ##############
### Array fastQC - Bowtie2 - sortBAM 
NUM=${#FILES[@]} # get size of array
ZBNUM=$(($NUM - 1))
if [ $ZBNUM -ge 0 ]; then 
  sbatch --wait --array=0-$ZBNUM%15 --output=${BaseFolder}ATACseq_%a.log --export=ALL,FILES="${BaseFolder}files.txt",BaseFolder=$BaseFolder,PairedEnd=$PairedEnd,Genome=$Genome  ${ScriptsFolder}ATAC_pipeline.sbs  
fi

wait