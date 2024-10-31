# R script to create count table

# args[1] = BaseFolder
# args[2] = PairedEnd
# args[3] = countMat
# args[4] = Genome
args=(commandArgs(TRUE))

print(args)

# Run

library(Rsubread)



BAMS<-dir(paste0(args[1], "BAM"),pattern=".sort.rmdup.rmblackls.rmchr.bam$", full.names=TRUE)
annot<-dir(paste0(args[1], "Analysis"), pattern = "Nahia_All.saf", full.names=TRUE)

# Single end
if(args[2]==0){
	Counts<-featureCounts(file=BAMS, annot.ext = annot, isGTFAnnotationFile = F,strandSpecific=0,nthreads=4,countMultiMappingReads=T)
}

# Paired end
if(args[2]==1){
	Counts<-featureCounts(file=BAMS, annot.ext = annot,isPairedEnd=TRUE)
}

save(Counts, file=args[3])
write.table(Counts$stat,file=paste0(args[1],"fastqc/FeatureCounts_stats_Nahia_All.summary"),row.names = F,quote = F,sep = "\t")
