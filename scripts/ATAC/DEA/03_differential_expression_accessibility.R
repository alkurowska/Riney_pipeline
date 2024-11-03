###########################################
####         ATAC-Seq pipeline         ####
####  03_differential_accessibility.R  ####
####              HUMAN                ####
###########################################

set.seed(1234)

# Libraries
library(edgeR)
library(limma)

###################
#### LOAD DATA ####
###################

# LOAD Cytogenetics
setwd("/ibex/user/kurowsaa/RINEY/FINAL/ATAC_data/Trying")
load("atac_data.RData")
dim(fish_atac) # 199 x 16

# remove HC young 
fish_atac <- fish_atac[fish_atac$Stage != "HC_young",] # 193 x 16

# Load ATAC data
setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/Count_matrix")
load("atac_to_dea.RData")
atac_data <- combat_edata

table(colnames(atac_data) == rownames(fish_atac))

# Relevant variables for the design matrix
stage = as.factor(fish_atac$Stage)
amp = as.factor(fish_atac$qall)
del1p = as.factor(fish_atac$del1p)
del17p = as.factor(fish_atac$del17p)
trans_14_16 = as.factor(fish_atac$trans_14_16)
trans_4_14 = as.factor(fish_atac$trans_4_14)

# design model
design <- model.matrix(~0 + stage + amp + del1p + del17p + trans_4_14 + trans_14_16)
colnames(design) <- gsub("stage", "", colnames(design))
colnames(design) <- gsub("amp", "", colnames(design))
colnames(design) <- gsub("del1p1p del", "del1p", colnames(design))
colnames(design) <- gsub("del17p17p del", "del17", colnames(design))
colnames(design)[8] <- "trans_4_14"
colnames(design)[9] <- "trans_14_16"


# Differential expression analysis
fit <- lmFit(atac_data, design)
fit2 <- eBayes(fit)

summary(decideTests(fit2))

# # Check the difference between each stage
myContrMatrix = makeContrasts(
   MM_SMM = MM - SMM, #trasitiion
   SMM_MGUS = SMM - MGUS, #trasitiion
   MM_MGUS = MM - MGUS, #trasitiion
   MGUS_HC = MGUS - HC, #trasitiion
   levels= colnames(design)
  )

## Get results
# estimate contrast for each gene
fit.contrasts <- contrasts.fit(fit, myContrMatrix)
fit.contrasts2 <- eBayes(fit.contrasts)

contrast_stats <- topTable(fit.contrasts2, number=nrow(atac_data), sort.by="none")

head(contrast_stats, 20)
length(which(contrast_stats$adj.P.Val < 0.05)) #19283
sum(contrast_stats$adj.P.Val < 0.01) #9109

dtest <- decideTests(fit.contrasts2)

contr.matrix <- myContrMatrix
all_padj <- matrix(ncol=ncol(contr.matrix), nrow=nrow(atac_data))
all_fch <- matrix(ncol=ncol(contr.matrix), nrow=nrow(atac_data))
sig <- matrix(0, ncol=ncol(contr.matrix), nrow=nrow(atac_data))

for(i in 1:ncol(contr.matrix)){
  target_stats <- topTable(fit.contrasts2, coef=i,number=nrow(atac_data), sort.by="none")
  
  all_padj[,i] <- target_stats[,"adj.P.Val"]
  all_fch[,i] <- target_stats[,"logFC"]
  sel <- intersect(which(all_padj[,i] < 0.05),which(all_fch[,i] > 0.58))
  sig[sel,i] <- 1
  sel2 <- intersect(which(all_padj[,i] < 0.05),which(all_fch[,i] < -0.58))
  sig[sel2,i] <- -1 
}

colnames(all_padj) = paste("p.adj_",colnames(contr.matrix),sep="")
colnames(all_fch) = paste("logFC_",colnames(contr.matrix),sep="")
colnames(sig) = paste("sig_",colnames(contr.matrix),sep="")

dea_results = data.frame(sig,all_padj,all_fch)
dim(dea_results)
rownames(dea_results) <- rownames(contrast_stats)
which(duplicated(rownames(dea_results)))

length(rownames(dea_results)[rowSums(abs(sig)!=0)>0]) #6971
list_genes <- rownames(dea_results)[apply(sig,1,FUN=function(x){return(sum(x!=0))})>0]
length(list_genes) #6971
#with a adj.pval<0.05 & |FC|>1.5 (logFC>0.58)
setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/DEA")
write.table(dea_results,"atac_dea_results.txt",sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)

# How many differential for each contrast
length(which(dea_results$sig_MM_SMM != 0)) # 0
length(which(dea_results$sig_SMM_MGUS != 0)) # 2
length(which(dea_results$sig_MGUS_HC != 0)) # 6362

table(dea_results$sig_MGUS_HC)
#    -1      0      1 
#  2204 133623   4158


table(dea_results$sig_SMM_MGUS)
# -1      0      1 
#  1 139983      1


# Plot logFC density

# density plot of logFC
toPlot <- data.frame(dea_results$logFC_MGUS_HC, dea_results$logFC_SMM_MGUS, dea_results$logFC_MM_SMM)
colnames(toPlot) <- c("T1", "T2", "T3")
rownames(toPlot) <- rownames(dea_results)

library(ggplot2)
png("densitylogFC.png", width=2000, height=2000, res=300)

ggplot(toPlot, aes(x = T1, fill = "T1")) +
    geom_density(alpha = 0.5) +
    geom_density(aes(x = T2, fill = "T2"), alpha = 0.5) +
    geom_density(aes(x = T3, fill = "T3"), alpha = 0.5) +
    xlab("log2FC") + 
    ylab("Density") +
    # dash line at zero
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    # add legend title
    scale_fill_manual(values = c("T1" = "red", "T2" = "green", "T3" = "blue"), labels = c("T1", "T2", "T3")) +
    guides(fill = guide_legend(title = "Transition", override.aes = aes(label = ""))) +
    labs(color = "") +
    ggtitle("Density of log2FC") +
    theme_minimal()
dev.off()


# Density plot of significant in T1
toPlot <- data.frame(dea_results$logFC_MGUS_HC[which(dea_results$sig_MGUS_HC != 0)], 
                     dea_results$logFC_SMM_MGUS[which(dea_results$sig_MGUS_HC != 0)], 
                     dea_results$logFC_MM_SMM[which(dea_results$sig_MGUS_HC != 0)])

colnames(toPlot) <- c("T1", "T2", "T3")
rownames(toPlot) <- rownames(dea_results)[which(dea_results$sig_MGUS_HC != 0)]

png("densitylogFC_sigT1.png", width=2000, height=2000, res=300)

ggplot(toPlot, aes(x = T1, fill = "T1")) +
    geom_density(alpha = 0.5) +
    geom_density(aes(x = T2, fill = "T2"), alpha = 0.5) +
    geom_density(aes(x = T3, fill = "T3"), alpha = 0.5) +
    xlab("log2FC") + 
    ylab("Density") +
    # dash line at zero
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    # add legend title
    scale_fill_manual(values = c("T1" = "red", "T2" = "green", "T3" = "blue"), labels = c("T1", "T2", "T3")) +
    guides(fill = guide_legend(title = "Transition", override.aes = aes(label = ""))) +
    labs(color = "") +
    ggtitle("Density of log2FC") +
    theme_minimal()
dev.off()


# Plot heatmap of significant peaks
toPlot <- unique(c(rownames(dea_results)[which(dea_results$sig_MGUS_HC != 0)],rownames(dea_results)[which(dea_results$sig_SMM_MGUS != 0)]))
count_table2 <- atac_data[toPlot,]
dim(count_table2)


#heatmap annotation
library(ComplexHeatmap)
library(RColorBrewer)

colors <- c(brewer.pal(n=8, name = 'Dark2'),brewer.pal(n=12, name = 'Set3'),brewer.pal(n=12, name = 'Paired'))
fish_anno <- fish_atac

fish_anno$translocation <- as.character(fish_anno$translocation)
fish_anno$translocation <- gsub(" ", "", fish_anno$translocation)

fish_anno$qall <- as.character(fish_anno$qall)
fish_anno[fish_anno$qall == '',]$qall <- "neutral"
fish_anno$del1p <- as.character(fish_anno$del1p)
fish_anno[fish_anno$del1p == '',]$del1p <- "neutral"
fish_anno$del17p <- as.character(fish_anno$del17p)
fish_anno[fish_anno$del17p == '',]$del17p <- "neutral"
fish_anno$qall <- gsub("qaber", "1q aberration", fish_anno$qall)

color_stage <- c("#33A02C", "#7570B3", "#E7298A", "#FF7F00") 
names(color_stage) <- c("HC", "MGUS", "MM", "SMM")

trans <- c("#8DD3C7", "#FB9A99","white")
names(trans) <- c("t(4;14)", "t(14;16)", "neutral")

del17 <- c("#FFED6F", "white")
names(del17) <- c("17p del", "neutral")

amp1q <- c("#A6CEE3", "white")
names(amp1q) <- c("1q aberration", "neutral")

del1p <- c("#666666", "white")
names(del1p) <- c("1p del", "neutral")

ha1 <- HeatmapAnnotation(
  Stage = fish_anno$Stage,
  Translocations = fish_anno$translocation,
  `17p deletion` = fish_anno$del17p,
  `1q aberration` = fish_anno$qall,
  `1p deletion` = fish_anno$del1p,
  col = list(Stage = as.factor(color_stage),
             Translocations = as.factor(trans),
             `17p deletion` = as.factor(del17),
             `1q aberration` = as.factor(amp1q),
             `1p deletion` = as.factor(del1p)),
  na_col = "#D9D9D9",
  show_annotation_name = TRUE)


data_to_plot2 <- t(scale(t(count_table2)))
dim(data_to_plot2)
getwd()


library(circlize)

pdf("atac_dea_heatmap.pdf", width=10, height=10)
Heatmap(data_to_plot2, name = "z-score", col = colorRamp2(c(-3, 0, 3), c("#4475B3", "white", "#D7342A")),
        top_annotation = ha1, show_row_dend = FALSE, show_column_dend = FALSE, 
        show_column_names = FALSE, show_row_names = FALSE, #row_names_gp = gpar(fontsize = 6), 
        cluster_rows = TRUE, column_split = factor(fish_anno$Stage, levels = c("HC", "MGUS", "SMM", "MM")),
        cluster_row_slices = FALSE, row_gap = unit(0, "mm"), row_title_rot = 0,
        cluster_column_slices = FALSE, column_gap = unit(0.5, "mm"),
        cluster_columns = TRUE, use_raster=TRUE)

dev.off()
