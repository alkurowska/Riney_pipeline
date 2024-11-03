# LOAD THE RAW DATA
setwd("/ibex/user/kurowsaa/RINEY/human/ATAC_data")
load("Rsubread_Counts_ATAC_Study1.Rdata") 
counts <- as.data.frame(Counts$counts) # 157817 x  248

# Load consensus peaks
setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/Count_matrix")
consensus_peaks <- read.table("consensus_peaks.txt", header = T, sep = "\t")

# check if the rownames of the counts matrix are the same as the consensus peaks
table(rownames(counts) == consensus_peaks$GeneID)
# TRUE
rownames(counts) <- consensus_peaks$name # change rownames of the count data

# Samples to use
setwd("/ibex/user/kurowsaa/RINEY/human_final/ATAC_data/Count_matrix")
samples_info <- read.table("atac_to_use_WyoungHC.txt", header = T)

counts <- counts[,colnames(counts) %in% samples_info$names] # 157817 x  199
dim(counts)

# change colnames of the count data
# Remove everything after "_S"
colnames(counts) <- gsub("_S.*", "", colnames(counts))

# Remove SEX chromosomes
counts <- counts[!grepl("chrX|chrY", rownames(counts)),] # 152912  x  199
dim(counts)

#####################
#### FRIP SCORES ####
####     BIN     ####
#####################

# Get the data from mapping
stats <- Counts$stat 
rownames(stats) <- stats$Status
stats <- stats[,-1] # remove first column with status names
colnames(stats) <- gsub("_S.*", "", colnames(stats))
stats <- as.data.frame(t(stats))
stats$Reads <- NA
stats$ConsFRIP <- NA

# Calculate the FRIP score of the consensus peaks
for (i in 1:nrow(stats)) {
  stats$Reads[i] <- sum(stats[i,1:14])
  stats$ConsFRIP[i] <- stats$Assigned[i]/stats$Reads[i]
}

stats$ConsFRIP <- round(stats$ConsFRIP,digits = 2)
# remove X from rownames
rownames(stats) <- gsub("X", "", rownames(stats))
stats$sample <- rownames(stats)

# Keep only the samples to use
stats <- stats[stats$sample %in% samples_info$sample,]

########
# Bins #
########

library(dplyr)

bins <- seq(0,1,0.05)

stats <- stats %>% mutate(bin = cut(ConsFRIP, breaks=bins,labels=F))
table(stats$bin)
# 2  3  4  5  6  7  8  9 10 11 14 
# 21 35 34 27 35 13 19  9  4  1  1

# All bins above 9 categorized as 9
for (i in 1:nrow(stats)) {
  if (isTRUE(stats$bin[i] > 9)) {
    stats$bin[i] <- 9
  }
}

# Assign bin group names
stats$bin <- factor(stats$bin,labels = c("2"="Group1","3"="Group2","4"="Group3","5"="Group4","6"="Group5","7"="Group6","8"="Group7","9"="Group8"))

rownames(samples_info) <- samples_info$sample
samples_info$bin <- stats[rownames(samples_info),]$bin


######################################
#### Counts distribution in peaks ####
####          PER STAGE           ####
######################################

# stage <- factor(samples_info$Stage, levels = c("HC_young", "HC", "MGUS", "SMM", "MM"))

# # Plot the distribution of the counts as histograms
# setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/Count_matrix")

# # Sum counts in each peak per stage
# toPlot <- data.frame(matrix(NA, nrow = nrow(counts), ncol = length(unique(stage))))
# colnames(toPlot) <- c("HC_young", "HC", "MGUS", "SMM", "MM")
# rownames(toPlot) <- rownames(counts)

# for (i in 1:ncol(toPlot)) {
#   # get sample names for the stage
#   samples <- samples_info[samples_info$Stage == colnames(toPlot)[i],]$sample

#   # sum counts in each peak for the samples in the stage
#   for (j in 1:nrow(toPlot)) {
#     toPlot[j,i] <- sum(counts[j,samples])
#   }
# }

# # Plot the distribution of the counts as histograms
# # x-axis is each peak in the order
# # y-axis is the sum of counts in each peak

# library(ggplot2)
# library(reshape2)

# toPlot$peak <- rownames(counts)
# toPlt <- melt(toPlot[1:100,])

# # Plot the distribution of the counts as line plots
# png("counts_distribution_per_stage_100.png", width = 3000, height = 600, res = 300)
# ggplot(data = toPlt[], aes(x = peak, y = value, color = variable, group = variable)) +
#   geom_line() +
#   labs(title = "Counts distribution per stage", x = "Peaks", y = "Counts") +
#   scale_color_manual(values = c("#ACD39E", "#33A02C", "#7570B3","#FF7F00", "#E7298A")) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
# dev.off()

#######################
#### Normalization ####
#######################

# Remove HC_young samples
counts <- counts[,samples_info[samples_info$Stage != "HC_young",]$sample]

summary(colnames(counts)==samples_info[samples_info$Stage != "HC_young",]$sample)
# TRUE
### create DGEList object
library(edgeR)
y <- DGEList(counts=counts)

## Filtering with the default parameters
stage <- factor(samples_info[samples_info$Stage != "HC_young",]$Stage, levels = c("HC", "MGUS", "SMM", "MM"))
keep <- filterByExpr(y, group=stage, min.count=10, min.total.count=15) 
table(keep)
# FALSE   TRUE 
# 12927 139985

# keep de peaks that pass the filter and recalculate the library size.
y <- y[keep, , keep.lib.sizes=FALSE]


## Normalization
# Data is normalized using the TMM function in EdgeR.
y2 <- calcNormFactors(y, method="TMM")
tmm_table <- cpm(y2)
setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/Count_matrix")
write.table(tmm_table, "tmm_data.txt", sep="\t", dec=".", quote=FALSE)


######## PCA ########
library(ggplot2)
library(RColorBrewer)

# Labels colors
color_stage <- c("#33A02C", "#7570B3","#FF7F00", "#E7298A") 
names(color_stage) <- c("HC", "MGUS", "SMM", "MM")

# PCA
# Raw counts
# log2 transformation
lcounts <- log2(y2$counts+1)
raw_pca <- prcomp(t(lcounts))

# TMM normalized counts
# log2 transformation
ltmm <- log2(tmm_table+1)
tmm_pca <- prcomp(t(ltmm))

library(factoextra)
png("PCA_variance_raw.png", width = 800, height = 600)
fviz_screeplot(raw_pca, addlabels = TRUE, ylim = c(0, 30))
dev.off()

png("PCA_variance_tmm.png", width = 800, height = 600)
fviz_screeplot(tmm_pca, addlabels = TRUE, ylim = c(0, 30))
dev.off()


# TMM PCA
var_prcomp <- tmm_pca$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

mds <- data.frame(PC1=tmm_pca$x[,1],PC2=tmm_pca$x[,2], Stage=stage)
labels <- rownames(mds)

p <- ggplot(mds, aes(x = PC1, y = PC2, color = Stage)) +
  scale_color_manual(values = color_stage) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  labs(x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 15, face = "bold"),
        legend.text=element_text(size=15))

ggsave(plot = p, dpi = 600, filename = "human_pca_tmm.jpg", width = 250, height = 200, units = "mm")


# RAW PCA
var_prcomp <- raw_pca$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

mds <- data.frame(PC1=raw_pca$x[,1],PC2=raw_pca$x[,2], Stage=stage)
labels <- rownames(mds)

p <- ggplot(mds, aes(x = PC1, y = PC2, color = Stage)) +
  scale_color_manual(values = color_stage) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  labs(x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 15, face = "bold"),
        legend.text=element_text(size=15))

ggsave(plot = p, dpi = 600, filename = "human_pca_raw.jpg", width = 250, height = 200, units = "mm")


# Color samples by FRIP score bins
# Labels colors
color_bin <- colorRampPalette(brewer.pal(8, "Set1"))(8)
names(color_bin) <- c("Group1", "Group2","Group3","Group4","Group5","Group6","Group7","Group8")


FRIP <- factor(samples_info[samples_info$Stage != "HC_young",]$bin, levels = c("Group1", "Group2","Group3","Group4","Group5","Group6","Group7","Group8"))

# TMM PCA
mds <- data.frame(PC1=tmm_pca$x[,1],PC2=tmm_pca$x[,2], FRIP=FRIP)
labels <- rownames(mds)

p <- ggplot(mds, aes(x = PC1, y = PC2, color = FRIP)) +
  scale_color_manual(values = color_bin) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  labs(x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 15, face = "bold"),
        legend.text=element_text(size=15))

ggsave(plot = p, dpi = 600, filename = "human_pca_tmm_bins.jpg", width = 250, height = 200, units = "mm")


# RAW PCA
mds <- data.frame(PC1=raw_pca$x[,1],PC2=raw_pca$x[,2], FRIP=FRIP)
labels <- rownames(mds)

p <- ggplot(mds, aes(x = PC1, y = PC2, color = FRIP)) +
  scale_color_manual(values = color_bin) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  labs(x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 15, face = "bold"),
        legend.text=element_text(size=15))

ggsave(plot = p, dpi = 600, filename = "human_pca_raw_bins.jpg", width = 250, height = 200, units = "mm")



######################
#### BATCH EFFECT ####
####     FRIP     ####
######################

library(sva)

# Create the model matrix

modcombat <- model.matrix(~Stage, data = samples_info[samples_info$Stage != "HC_young",])
batch <- samples_info[samples_info$Stage != "HC_young",]$bin

#input_combat <- norm_data
input_combat <- log2(tmm_table+1)
combat_edata <- ComBat(dat = input_combat, batch = batch, mod = modcombat, par.prior = T, mean.only = T)


# PCA
# Combat
combat_pca <- prcomp(t(combat_edata))


png("PCA_variance_combat.png", width = 800, height = 600)
fviz_screeplot(combat_pca, addlabels = TRUE, ylim = c(0, 30))
dev.off()

# Combat PCA
var_prcomp <- combat_pca$sdev^2
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

mds <- data.frame(PC1=combat_pca$x[,1],PC2=combat_pca$x[,2], Stage=stage)
labels <- rownames(mds)

p <- ggplot(mds, aes(x = PC1, y = PC2, color = Stage)) +
  scale_color_manual(values = color_stage) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  labs(x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 15, face = "bold"),
        legend.text=element_text(size=15))

ggsave(plot = p, dpi = 600, filename = "human_pca_combat.jpg", width = 250, height = 200, units = "mm")


# Color samples by FRIP score bins

mds <- data.frame(PC1=combat_pca$x[,1],PC2=combat_pca$x[,2], FRIP=FRIP)
labels <- rownames(mds)

p <- ggplot(mds, aes(x = PC1, y = PC2, color = FRIP)) +
  scale_color_manual(values = color_bin) + 
  geom_point(size = 6, show.legend = TRUE, shape = 19) +
  labs(x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 20, face = "bold"),legend.title = element_text(size = 15, face = "bold"),
        legend.text=element_text(size=15))

ggsave(plot = p, dpi = 600, filename = "human_pca_combat_bins.jpg", width = 250, height = 200, units = "mm")


# Save data
metadata <- samples_info[samples_info$Stage != "HC_young",]
setwd("/ibex/user/kurowsaa/RINEY/NEW/ATAC/Count_matrix")
save(combat_edata, metadata, file = "atac_to_dea.RData")


