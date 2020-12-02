# Written by Emily Hughes, May 7 2020
# This code was written as part of the Boston University Bioinformatics Masters Program
# BF528: Applications in Translational Bioinformatics, Project 5

# This code was used to replicate the biologist role of Project 2- replicating O'Meara et al
library(ggplot2)
library(dplyr)
library(ggpubr)

# Set the working directory
#setwd("~/Documents/Applications 2/Project_5/proj2_biologist")

# Part 7-1: Compare trends of genes in FPKM table with those of Figure 1D in paper
# Note, sample data will be used to complete this analysis

# Read in the data
fpkm_matrix <- read.table('fpkm_matrix.csv', header = TRUE)

# Subset the matrix
# Get averages for each value
fpkm_matrix$Ad <- (fpkm_matrix$Ad_1_FPKM + fpkm_matrix$Ad_2_FPKM)/2
fpkm_matrix$P4 <- (fpkm_matrix$P4_1_FPKM + fpkm_matrix$P4_2_FPKM)/2
fpkm_matrix$P7 <- (fpkm_matrix$P7_1_FPKM + fpkm_matrix$P7_2_FPKM)/2
# Get only these averages
fpkm_matrix <- fpkm_matrix[ , !(names(fpkm_matrix) %in% c("Ad_1_FPKM", "Ad_2_FPKM", "P4_1_FPKM", "P4_2_FPKM", "P7_1_FPKM", "P7_2_FPKM"))]

# Mitochondria data
genes <- list('ENSMUSG00000023861'='Mpc1','ENSMUSG00000024997'='Prdx3', 'ENSMUSG00000032047' ='Acat1',  'ENSMUSG00000025465'='Echs1',  'ENSMUSG00000014606'='Slc25a11',  
'ENSMUSG00000026664'='Phyh')
mito <- subset(fpkm_matrix, fpkm_matrix$tracking_id %in% names(genes))
for(id in names(genes)){
  for(row in 1:6){
    if(mito$tracking_id[row] == id){
      mito$gene_name[row] <- genes[id]
    }
  }
}

# Clean up
names(mito)[2] <- "P0"
mito <- mito[, 2:6]
clean <- mito %>% gather(time, expression, c("P0", "P4", "P7", "Ad"), -gene_name)
df2 <- as.data.frame(lapply(clean, unlist))
df2 <- df2[order(df2[, 1]), ]
df2$time <- factor(df2$time, levels=c("P0", "P4", "P7", "Ad"))

a <- ggplot(df2, aes(x=time, y=expression, group=gene_name)) + geom_point(aes(color=gene_name)) + 
  geom_line(aes(color=gene_name)) + labs(title="Mitochondria",x="", y = "FPKM") +
  theme(legend.title = element_blank())

# Sarcomere data
genes <- list("ENSMUSG00000028273"='Pdlim5','ENSMUSG00000032648'='Pygm', 'ENSMUSG00000028116' ='Myoz2',  'ENSMUSG00000026208'='Des',  'ENSMUSG00000030470'='Csrp3',  
              'ENSMUSG00000007877'='Tcap', 'ENSMUSG00000032060'='Cryab')
sarc <- subset(fpkm_matrix, fpkm_matrix$tracking_id %in% names(genes))
sarc[,'gene_name'] <- NA
for(id in names(genes)){
  for(row in 1:7){
    if(sarc$tracking_id[row] == id){
      sarc$gene_name[row] <- genes[id]
    }
  }
}
# Clean up
names(sarc)[2] <- "P0"
sarc <- sarc[, 2:6]
clean <- sarc %>% gather(time, expression, c("P0", "P4", "P7", "Ad"), -gene_name)
df2 <- as.data.frame(lapply(clean, unlist))
df2 <- df2[order(df2[, 1]), ]
df2$time <- factor(df2$time, levels=c("P0", "P4", "P7", "Ad"))

b <- ggplot(df2, aes(x=time, y=expression, group=gene_name)) + geom_point(aes(color=gene_name)) + 
  geom_line(aes(color=gene_name)) + labs(title="Sarcomere",x="", y = "FPKM") +
  theme(legend.title = element_blank())

# Cell cycle
genes <- list("ENSMUSG00000029283"='Cdc7','ENSMUSG00000046179'='E2f8', 'ENSMUSG00000069089' ='Cdk7',  'ENSMUSG00000066149'='Cdc26',  'ENSMUSG00000017499'='Cdc6',  
              'ENSMUSG00000027490'='E2f1', 'ENSMUSG00000020687'='Cdc27', 'ENSMUSG00000000028' = 'Cdc45', 'ENSMUSG00000027323' = 'Rad51', 'ENSMUSG00000020897'='Aurkb',
              'ENSMUSG00000024370'='Cdc23', 'ENSMUSG00000022070'='Bora')
cell <- subset(fpkm_matrix, fpkm_matrix$tracking_id %in% names(genes))
cell[,'gene_name'] <- NA
for(id in names(genes)){
  for(row in 1:12){
    if(cell$tracking_id[row] == id){
      cell$gene_name[row] <- genes[id]
    }
  }
}
# Clean up
names(cell)[2] <- "P0"
cell <- cell[, 2:6]
clean <- cell %>% gather(time, expression, c("P0", "P4", "P7", "Ad"), -gene_name)
df2 <- as.data.frame(lapply(clean, unlist))
df2 <- df2[order(df2[, 1]), ]
df2$time <- factor(df2$time, levels=c("P0", "P4", "P7", "Ad"))

c <- ggplot(df2, aes(x=time, y=expression, group=gene_name)) + geom_point(aes(color=gene_name)) + 
  geom_line(aes(color=gene_name)) + labs(title="Cell cycle",x="", y = "FPKM") +
  theme(legend.title = element_blank())

ggarrange(b, a, c, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

# Part 7-3: Creating a heat map
ad_1 <- read.table('Ad_1.fpkm_tracking', header = TRUE)
ad_2 <- read.table('Ad_2.fpkm_tracking', header = TRUE)
p0_2 <- read.table('P0_2.fpkm_tracking', header = TRUE)
p4_1 <- read.table('P4_1.fpkm_tracking', header = TRUE)
p4_2 <- read.table('P4_2.fpkm_tracking', header = TRUE)
p7_1 <- read.table('P7_1.fpkm_tracking', header = TRUE)
p7_2 <- read.table('P7_2.fpkm_tracking', header = TRUE)
diff <- read.table('gene_exp.diff', header=TRUE)

# Concatenate the FPKM columns of each of the tracking datasets into a single dataframe
track_id <- as.data.frame(p0_2$tracking_id)
gene <- as.data.frame(p0_2$gene_short_name)
full <- cbind(track_id, gene)
colnames(full) <- c('tracking_id', 'gene_name')

add_cols <- function(df, big, name){
  new <- list()
  for(row in 1:37468){
    same <- match(big$tracking_id[row], df$tracking_id)
    new <- c(new, df$FPKM[same])
  }
  new <- t(as.data.frame(new))
  final <- cbind(big, new)
  colnames(final)[-1] <- name
  return(final)
}

full <- add_cols(p0_2, full, 'P0')
full <- add_cols(p4_1, full, 'P4_1')
full <- add_cols(p4_2, full, 'P4_2')
full <- add_cols(p7_1, full, 'P7_1')
full <- add_cols(p7_2, full, 'P7_2')
full <- add_cols(ad_1, full, 'Ad_1')
full <- add_cols(ad_2, full, 'Ad_2')
colnames(full) <- c('tracking_id', 'gene_name', 'P0_2', "P4_1", 'P4_2', 'P7_1', 'P7_2', "Ad_1", "Ad_2")

# Get top 500 differentially expressed genes between P4 and P7 (example data)
# get rid of infinite values
diff <- subset(diff, (diff$log2.fold_change. != "Inf" & diff$log2.fold_change. != "-Inf"))
genes <- top_n(diff, 500, abs(log2.fold_change.)) # get top 1000 genes of highest absolute log fold change
gene_names <- genes$gene

clean <- subset(full, full$gene_name %in% gene_names)
clean <- as.matrix(clean)
really_clean <- clean[, 3:9]
really_clean <- really_clean[apply(really_clean[,-1], 1, function(x) !all(x<0.0000001)),]
mode(really_clean) <- 'numeric'
rownames(really_clean) <- c()
heatmap(really_clean, scale='row')