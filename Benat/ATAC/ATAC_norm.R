
#### ATAC GLOBAL ANALYSIS ####

ATAC_Counts <- Counts$counts
ATAC_Counts <- as.data.frame(ATAC_Counts)
colnames(ATAC_Counts)

ATAC_stats <- Counts$stat
row.names(ATAC_stats) <- ATAC_stats$Status
ATAC_stats <- ATAC_stats[,-1]
ATAC_stats <- t(ATAC_stats)
ATAC_stats <- as.data.frame(ATAC_stats)
ATAC_stats$Reads <- NA
ATAC_stats$ConsFRIP <- NA
for (i in 1:nrow(ATAC_stats)) {
  ATAC_stats$Reads[i] <- sum(ATAC_stats[i,1:14])
  ATAC_stats$ConsFRIP[i] <- ATAC_stats$Assigned[i]/ATAC_stats$Reads[i]
}

ATAC_stats$ConsFRIP <- round(ATAC_stats$ConsFRIP,digits = 2)
ATAC_stats$Sample <- row.names(ATAC_stats)

for (i in 1:nrow(ATAC_stats)) {
  ATAC_stats$Sample[i] <- unlist(strsplit(row.names(ATAC_stats)[i],split = "_S"))[1]
}

samples_info$ConsFRIP <- NA
for (i in 1:nrow(ATAC_stats)) {
  for (j in 1:nrow(samples_info)) {
    if (isTRUE(ATAC_stats$Sample[i] == samples_info$Sample[j])) {
      samples_info$ConsFRIP[j] <- ATAC_stats$ConsFRIP[i]
    }
  }
}

samples_info_ST1 <- samples_info[samples_info$Disease %in% c("MM","MM-HR","healthy donors young","AL amyloidosis"),]
samples_info_ST1$Disease <- gsub("MM-HR","MM",samples_info_ST1$Disease)

library(dplyr)

bins <- seq(0,1,0.05)

samples_info_ST1 <- samples_info_ST1 %>% mutate(bin = cut(ConsFRIP, breaks=bins,labels=F))
table(samples_info_ST1$Disease,samples_info_ST1$bin)

for (i in 1:nrow(samples_info_ST1)) {
  if (isTRUE(samples_info_ST1$bin[i] > 9)) {
    samples_info_ST1$bin[i] <- 9
  }
}

samples_info_ST1$bin <- factor(samples_info_ST1$bin,labels = c("2"="Group1","3"="Group2","4"="Group3","5"="Group4","6"="Group5","7"="Group6","8"="Group7","9"="Group8"))


table(samples_info_ST1$Disease,samples_info_ST1$bin)
unique(samples_info_ST1$bin)

for (i in 1:ncol(ATAC_Counts)) {
  colnames(ATAC_Counts)[i] <- unlist(strsplit(colnames(ATAC_Counts)[i],split = "_S"))[1]
}

ATAC_Counts <- ATAC_Counts[,colnames(ATAC_Counts) %in% samples_info_ST1$Sample]

orden <- samples_info_ST1$Sample
ATAC_Counts <- ATAC_Counts[,orden]

all(colnames(ATAC_Counts) == samples_info_ST1$Sample)

which(colnames(ATAC_Counts)=="AC-BM616")
ATAC_Counts <- ATAC_Counts[,-which(colnames(ATAC_Counts)=="AC-BM616")]
samples_info_ST1 <- samples_info_ST1[-which(colnames(ATAC_Counts)=="AC-BM616"),]

Sex_CHR <- which(Counts$annotation$Chr %in% c("chrX","chrY"))
ATAC_Counts <- ATAC_Counts[-c(Sex_CHR),]

colnames_interes <- colnames(BBDD_Final_Componentes)[9:22]

new_anno <- as.data.frame(matrix(NA, nrow = nrow(samples_info_ST1), ncol = length(colnames_interes)))
colnames(new_anno) <- colnames_interes
samples_info_ST1 <- cbind(samples_info_ST1,new_anno)

for (i in 1:nrow(samples_info_ST1)) {
  for (j in 1:nrow(BBDD_Final_Componentes)) {
    if (isTRUE(samples_info_ST1$Sample[i] == BBDD_Final_Componentes$Sample[j])) {
      samples_info_ST1[i,9:22] <- BBDD_Final_Componentes[j,9:22]
    }
  }
}

which(is.na(samples_info_ST1$Status) & samples_info_ST1$Disease == "MM")



all(colnames(ATAC_Counts) == samples_info_ST1$Sample)

y_ATAC <- DGEList(ATAC_Counts,samples = samples_info_ST1,group = samples_info_ST1$Disease,genes = Counts$annotation[-(Sex_CHR),])
y_ATAC <- calcNormFactors(y_ATAC,method = "TMM")
keep <- filterByExpr(y_ATAC,group = samples_info_ST1$Disease)

y_ATAC <- y_ATAC[keep,,keep.lib.sizes=T]

modcombat <- model.matrix(~Disease, data = samples_info_ST1)
batch <- samples_info_ST1$bin

library(sva)

input_combat <- log2(y_ATAC$counts+1)
combat_edata <- ComBat(dat = input_combat, batch = batch, mod = modcombat, par.prior = T, mean.only = T)
