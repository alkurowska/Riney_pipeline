# Consensus peak
# 
library(rtracklayer)
library(GenomicRanges)
library(ChIPpeakAnno)
library(limma)
library(readxl)

`%ni%` <- Negate(`%in%`)

peak_files<-list.files(path = "/Users/bariceta/hpc/RINEY/ATAC/OLD/PeakCalling", 
                       pattern = ".*_peaks.narrowPeak", full.names = T)
peak_files<-c(peak_files, list.files(path = "/Users/bariceta/hpc/RINEY/ATAC/PeakCalling", 
                                     pattern = ".*_peaks.narrowPeak", full.names = T))

aux <- rep(NA,248)

for (i in 1:248) {
  tmp <- unlist(strsplit(peak_files[i],split="/"))[8]
  tmp <- unlist(strsplit(tmp,split="_S"))[1]
  aux[i] <- tmp
}
for (i in 206:219) {
  tmp <- unlist(strsplit(peak_files[i],split="/"))[5]
  tmp <- unlist(strsplit(tmp,split="_S"))[1]
  aux[i] <- tmp
}

BBDD <- read_excel("BBDD.xlsx", sheet = "bulk")

ATAC_BBDD <- as.data.frame(matrix(NA,nrow = 248,ncol = 2))
colnames(ATAC_BBDD) <- c("Sample","Disease")

ATAC_BBDD$Sample <- aux

for (i in 1:248) {
  for (j in 1:297) {
    if (isTRUE(ATAC_BBDD$Sample[i] == BBDD$Sample[j])) {
      ATAC_BBDD$Disease[i] <- BBDD$Disease[j]
    }
  }
}

ATAC_BBDD$Disease[3] <- "SMM"
ATAC_BBDD$Disease[29] <- "MM"
ATAC_BBDD$Disease[30] <- "MM"
ATAC_BBDD$Disease[56] <- "SMM-HR"
ATAC_BBDD$Disease[59] <- "SMM-HR"
ATAC_BBDD$Disease[83] <- "SMM-HR"
ATAC_BBDD$Disease[84] <- "SMM"



ATAC_BBDD$Disease[97] <- "SMM"


ATAC_BBDD$Disease[185] <- "AL amyloidosis"
ATAC_BBDD$Disease[186] <- "AL amyloidosis"
ATAC_BBDD$Disease[187] <- "MM&AL"
ATAC_BBDD$Disease[189] <- "AL amyloidosis"
ATAC_BBDD$Disease[192] <- "AL amyloidosis"
ATAC_BBDD$Disease[193] <- "AL amyloidosis"
ATAC_BBDD$Disease[194] <- "AL amyloidosis"
ATAC_BBDD$Disease[196] <- "AL amyloidosis"
ATAC_BBDD$Disease[202] <- "BICLONAL&AL"
ATAC_BBDD$Disease[203] <- "AL amyloidosis"
ATAC_BBDD$Disease[207] <- "AL amyloidosis"
ATAC_BBDD$Disease[210] <- "AL amyloidosis"
ATAC_BBDD$Disease[211] <- "MM&AL"
ATAC_BBDD$Disease[233] <- "healthy donors young"
ATAC_BBDD$Disease[235] <- "healthy donors young"


ATAC_BBDD[-131,]


MM <- which(ATAC_BBDD$Disease %in% c("MM","MM-HR"))
SMM <- which(ATAC_BBDD$Disease %in% c("SMM","SMM-HR"))
MGUS <- which(ATAC_BBDD$Disease == "MGUS")
AL <- which(ATAC_BBDD$Disease == "AL amyloidosis")
HealthyYoung <- which(ATAC_BBDD$Disease == "healthy donors young")
HealthyOld <- which(ATAC_BBDD$Disease == "healthy donors")
Misc <- which(ATAC_BBDD$Disease %ni% c("MM","MM-HR","SMM","SMM-HR","MGUS","AL amyloidosis","healthy donors young"))

# peak files per condition

peak_files_MM<-peak_files[MM]
peak_files_SMM<-peak_files[SMM]
peak_files_MGUS<-peak_files[MGUS]
peak_files_AL<-peak_files[AL]
peak_files_HDY<-peak_files[HealthyYoung]
peak_files_HDO<-peak_files[HealthyOld]


# Consensus MM #

extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

peak_MM_granges<-lapply(peak_files_MM, import, extraCols = extraCols_narrowPeak, format = "BED")
summary(peak_MM_granges)

peak_MM_grangeslist<-GRangesList(peak_MM_granges)
peak_MM_grangeslist


# Select regions covered by at least 2 samples and merge peaks closer than 30bp

peak_coverage_MM<-coverage(peak_MM_grangeslist)
peak_coverage_MM
covered_ranges_MM<-slice(x=peak_coverage_MM, lower =2, rangesOnly = T)
covered_ranges_MM
covered_granges_MM<-GRanges(covered_ranges_MM)
summary(covered_granges_MM)
reduced_covered_granges_MM<-reduce(covered_granges_MM, min.gapwidth=30)
summary(reduced_covered_granges_MM)
export(object = reduced_covered_granges_MM, "consensus_MM_2_Final.bed")



# Consensus SMM #


extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

peak_SMM_granges<-lapply(peak_files_SMM, import, extraCols = extraCols_narrowPeak, format = "BED")
summary(peak_SMM_granges)

peak_SMM_grangeslist<-GRangesList(peak_SMM_granges)
peak_SMM_grangeslist


# Select regions covered by at least 2 samples and merge peaks closer than 30bp

peak_coverage_SMM<-coverage(peak_SMM_grangeslist)
peak_coverage_SMM
covered_ranges_SMM<-slice(x = peak_coverage_SMM, lower =2, rangesOnly = T)
covered_ranges_SMM
covered_granges_SMM<-GRanges(covered_ranges_SMM)
summary(covered_granges_SMM)
reduced_covered_granges_SMM<-reduce(covered_granges_SMM, min.gapwidth=30)
summary(reduced_covered_granges_SMM)
export(object = reduced_covered_granges_SMM, "consensus_SMM_2_Final.bed")



# Consensus MGUS #

extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

peak_MGUS_granges<-lapply(peak_files_MGUS, import, extraCols = extraCols_narrowPeak, format = "BED")
summary(peak_MGUS_granges)

peak_MGUS_grangeslist<-GRangesList(peak_MGUS_granges)
peak_MGUS_grangeslist


# Select regions covered by at least 2 samples and merge peaks closer than 30bp

peak_coverage_MGUS<-coverage(peak_MGUS_grangeslist)
peak_coverage_MGUS
covered_ranges_MGUS<-slice(x = peak_coverage_MGUS, lower = 2, rangesOnly = T)
covered_ranges_MGUS
covered_granges_MGUS<-GRanges(covered_ranges_MGUS)
summary(covered_granges_MGUS)
reduced_covered_granges_MGUS<-reduce(covered_granges_MGUS, min.gapwidth=30)
summary(reduced_covered_granges_MGUS)
export(object = reduced_covered_granges_MGUS, "consensus_MGUS_2_Final.bed")


# Consensus AL #

extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

peak_AL_granges<-lapply(peak_files_AL, import, extraCols = extraCols_narrowPeak, format = "BED")
summary(peak_AL_granges)

peak_AL_grangeslist<-GRangesList(peak_AL_granges)
peak_AL_grangeslist


# Select regions covered by at least 2 samples and merge peaks closer than 30bp

peak_coverage_AL<-coverage(peak_AL_grangeslist)
peak_coverage_AL
covered_ranges_AL<-slice(x = peak_coverage_AL, lower = 2, rangesOnly = T)
covered_ranges_AL
covered_granges_AL<-GRanges(covered_ranges_AL)
summary(covered_granges_AL)
reduced_covered_granges_AL<-reduce(covered_granges_AL, min.gapwidth=30)
summary(reduced_covered_granges_AL)
export(object = reduced_covered_granges_AL, "consensus_AL_2_Final.bed")


# Consensus Healthy #

extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

peak_HDY_granges<-lapply(peak_files_HDY, import, extraCols = extraCols_narrowPeak, format = "BED")
summary(peak_HDY_granges)

peak_HDY_grangeslist<-GRangesList(peak_HDY_granges)
peak_HDY_grangeslist


# Select regions covered by at least 1 sample and merge peaks closer than 30bp

peak_coverage_HDY<-coverage(peak_HDY_grangeslist)
peak_coverage_HDY
covered_ranges_HDY<-slice(peak_coverage_HDY, lower = 2, rangesOnly = T)
covered_ranges_HDY
covered_granges_HDY<-GRanges(covered_ranges_HDY)
summary(covered_granges_HDY)
reduced_covered_granges_HDY<-reduce(covered_granges_HDY, min.gapwidth=30)
summary(reduced_covered_granges_HDY)
export(object = reduced_covered_granges_HDY, "consensus_HDY_1_Final.bed")


# Consensus Healthy OLD#
peak_files_HDO <- peak_files_HDO[c(1:2,4:9)]

extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

peak_HDO_granges<-lapply(peak_files_HDO, import, extraCols = extraCols_narrowPeak, format = "BED")
summary(peak_HDO_granges)

peak_HDO_grangeslist<-GRangesList(peak_HDO_granges)
peak_HDO_grangeslist


# Select regions covered by at least 2 sample and merge peaks closer than 30bp

peak_coverage_HDO<-coverage(peak_HDO_grangeslist)
peak_coverage_HDO
covered_ranges_HDO<-slice(peak_coverage_HDO, lower = 2, rangesOnly = T)
covered_ranges_HDO
covered_granges_HDO<-GRanges(covered_ranges_HDO)
summary(covered_granges_HDO)
reduced_covered_granges_HDO<-reduce(covered_granges_HDO, min.gapwidth=30)
summary(reduced_covered_granges_HDO)
export(object = reduced_covered_granges_HDO, "consensus_HDO_2_Final.bed")

# FINAL CONSENSUS: Select all peaks (common and unique) and merge peaks closer than Xbp

peak_overlap<-findOverlapsOfPeaks(reduced_covered_granges_MM,reduced_covered_granges_SMM, reduced_covered_granges_MGUS,
                                  reduced_covered_granges_AL,reduced_covered_granges_HDO,
                                  maxgap = 30, connectedPeaks = "merge")

peak_overlap3 <- findOverlapsOfPeaks(peak_overlap,reduced_covered_granges_HDY,maxgap = 30,connectedPeaks = "merge")
(peak_overlap,reduced_covered_granges_HDY,maxgap = 30,connectedPeaks = "merge")

peak_overlap <- findOverlapsOfPeaks(reduced_covered_granges_MM,reduced_covered_granges_HDO,
                                    maxgap = 30, connectedPeaks = "merge")

vennDiagram(peak_overlap$venn_cnt, names = c("MM","SMM", "MGUS","AL","Elderly Donors"),circle.col = c("blue","green","pink","red","purple"))

vennDiagram(peak_overlap$venn_cnt,names = c("MM","Elderly Donors"),circle.col = c("blue","green"))
peak_overlap3<-peak_overlap3$peaklist
peak_overlap<-c(peak_overlap[[1]], peak_overlap[[2]], peak_overlap[[3]], peak_overlap[[4]],
                peak_overlap[[5]], peak_overlap[[6]], peak_overlap[[7]],peak_overlap[[8]],
                peak_overlap[[9]], peak_overlap[[10]], peak_overlap[[11]], peak_overlap[[12]],
                peak_overlap[[13]], peak_overlap[[14]], peak_overlap[[15]],peak_overlap[[16]],
                peak_overlap[[17]], peak_overlap[[18]], peak_overlap[[19]], peak_overlap[[20]],
                peak_overlap[[21]], peak_overlap[[22]], peak_overlap[[23]],peak_overlap[[24]],
                peak_overlap[[25]],peak_overlap[[26]],peak_overlap[[27]])
peak_overlap3 <- c(peak_overlap3[[1]],peak_overlap3[[2]],peak_overlap3[[3]])
peak_overlap<-sort(peak_overlap3)
summary(peak_overlap)
head(peak_overlap)
consensus<-as.data.frame(peak_overlap)
consensus<-consensus[,1:3]
consensus[,4]<-paste0("ConsensusPeak_", 1:nrow(consensus))
consensus[,5]<-0
consensus[,6]<-"."
write.table(x = consensus, file = "consensus_2_new.bed", quote = F, col.names = F, row.names = F, sep = "\t")


# -------------------------------------------------

# convert to SAF format - GeneID, Chr, Start, End and Strand
consensus<-read.table("consensus_2_new.bed")
rownames(consensus)<-paste0("ConsensusPeak_", 1:nrow(consensus))
consensus<-consensus[,1:4]
consensus$Strand<-"*"
colnames(consensus)<-c("Chr", "Start", "End", "GeneID", "Strand")
save(consensus, file = "consensus_2.Rdata")

saf <- consensus[,c(4,1:3,5)]
write.table(x=saf,file = "Consensus_ATAC.saf",quote = F,col.names = T,row.names = F,sep = "\t")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

annotated <- annotatePeak(peak_overlap,TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)

ChIPseeker::plotAnnoBar(annotated)
