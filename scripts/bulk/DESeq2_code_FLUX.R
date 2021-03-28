source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("apeglm")
biocLite("vsn")
biocLite("pcaExplorer")


library("png")
library("QoRTs")
library("JunctionSeq")
library("DESeq2")
library("DEXSeq")
library("edgeR")
library("data.table")
library("qdapTools")
library("pheatmap")
library("tidyverse")

#set working directory
setwd("/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/analysis_r/DESeq2/")

#set sample directory names
run2245_dir_prefix = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/run2245_output/featureCounts"
run2285_dir_prefix = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/run2285_output/featureCounts"

#testing that geneid list in the same between samples
test <- read.table(file = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/run2245_output/featureCounts/Sample_103916/Sample_103916+=_008_featureCounts.txt", 
                   sep = '\t',
                   header = TRUE,
                   stringsAsFactors = FALSE)
test2 <- read.table(file = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/run2285_output/featureCounts/Sample_105820/Sample_105820+=_featureCounts.txt", 
                    sep = '\t',
                    header = TRUE,
                    stringsAsFactors = FALSE)

#geneid tests TRUE for between 2245/2285 and different samples within 2245
identical(test$Geneid, test2$Geneid)


#############################################
# Building featureCounts expression matrix
#############################################

#MUST RUN THIS WHOLE BLOCK IN SEQUENCE
#intialize name and count vectors the expression count matrix
Geneid <- test$Geneid
counts <- data.frame(Geneid, stringsAsFactors = FALSE)
decoder_i <- 1

#cbind run2245 count information

for (run2245_id in list.files(path = run2245_dir_prefix,
                                 pattern = ".txt$",
                                 recursive = TRUE,
                                 full.names = TRUE)) {
  
  #read featureCount file into a dataframe
  featureCount.mtx <- read.table(run2245_id, sep = '\t', header = T, stringsAsFactors = FALSE)
  
  #retrieve the sample id using placenta_sort_decoder, developed for QoRTs - lists sample information
  #beginning with the first run in 2245 through the end of run2285; MUST RESET decoder_i before running loop
  sample_name <- placenta_sort_decoder_df$unique.ID[decoder_i]
  decoder_i <- decoder_i + 1
  
  #subset the dataframe to just the expression counts
  featureCount.mtx <- featureCount.mtx[,7]
  
  #rename the column with unique.ID
  names(featureCount.mtx)[names(featureCount.mtx) == featureCount.mtx] <- sample_name
  
  #cbind the the list of geneIDs to the expression count matrix of the sample
  counts <- cbind.data.frame(counts, featureCount.mtx)
}

#append run2285 information to vectors
for (run2285_id in list.files(path = run2285_dir_prefix,
                              pattern = ".txt$",
                              recursive = TRUE,
                              full.names = TRUE)) {
  
  #read featureCount file into a dataframe
  featureCount.mtx <- read.table(run2285_id, sep = '\t', header = T, stringsAsFactors = FALSE)
  
  #retrieve the sample id using placenta_sort_decoder, developed for QoRTs - lists sample information
  #beginning with the first run in 2245 through the end of run2285; MUST RESET decoder_i before running loop
  #this loop WILL NOT work without running the run2245 appending loop above
  sample_name <- placenta_sort_decoder_df$unique.ID[decoder_i]
  decoder_i <- decoder_i + 1
  
  #subset the dataframe to just the expression counts
  featureCount.mtx <- featureCount.mtx[,7]
  
  #rename the column with unique.ID
  names(featureCount.mtx)[names(featureCount.mtx) == featureCount.mtx] <- sample_name
  
  #cbind the the list of geneIDs to the expression count matrix of the sample
  counts <- cbind.data.frame(counts, featureCount.mtx)
}

#rename each column because it didn't work above
for(i in (2:length(colnames(counts)))) {
  colnames(counts)[i] <- placenta_sort_decoder_df$unique.ID[i-1]
}


#write.table( cts, file = "placenta_sort_counts.txt", sep = "\t", row.names = FALSE)

###################################################
# Entering the data into DESeq2 and Vignette
###################################################
cts <- read.table(file = "placenta_sort_counts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#these commands were used to take the output from line 101 and push column 1 into rownames (gene IDs)
#rownames(cts) <- cts[,1]
#cts <- cts[,-1]

coldata = read.table(file = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/analysis_r/placenta_sort_decoder.txt",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = FALSE)

#modifying qc.data.prefix to exclude lane identifer for collapseReplicates in DESeq2
for( j in 1:length(coldata$qc.data.prefix)) {
  coldata$qc.data.prefix[j] <- substr(coldata$qc.data.prefix[j] , 1, 13)
}
colnames(coldata)[8] = "run.sample.ID"

#save as a new file
write.table( coldata, file = "col_data_for_DESeq2.txt", sep = "\t", row.names = FALSE)

#read in new coldata for DESeq2
coldata = read.table(file = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/analysis_r/DESeq2/col_data_for_DESeq2.txt",
                     sep = "\t",
                     header = TRUE,
                     stringsAsFactors = TRUE)

#lane.ID and/or is.Paired were linear combinations of eacher and so couldn't run in the model
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ sample.ID + is.Paired + cell.type)

#################################################################
# Collapse Technical Replicates
#################################################################

dds <- collapseReplicates(dds, dds$run.sample.ID, renameCols = TRUE)

#################################################################
# Minimal Pre-filtering
#################################################################

#keep only rows that have more than 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

###############################################################################
# More Pre-filtering based on https://support.bioconductor.org/p/65091/ thread
###############################################################################

#implemented due to convergence problems in model fitting
#require normalized count of at least 10 in three samples
dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
dds <- dds[filter,]


#################################################################
# Setting factor levels (alphabetical by default)
#################################################################
dds$sample.ID <- factor(dds$sample.ID, levels = c("KC36", "KC39", "KC40", "KC48"))
dds$is.Paired <- factor(dds$is.Paired, levels = c("single", "paired"))
#composite is the reference group for cell type as a result

#################################################################
# Run DEX Analysis
#################################################################
dds <- DESeq(dds)
res <- results(dds)
res

#log output error for uncollapsed cell type model:
#fitting model and testing
#58 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest

#log output error for collapsed cell type model:
#fitting model and testing
#47 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest

#################################################################
# Specify a contrast between cell types
#################################################################

res <- results(dds, constrast = c("composite", "syncytiotrophoblast, leukocyte",
                                  "hofbauer", "fibroblast", "cytotrophoblast",
                                  "leukocyte", "extravillous trophoblast"))
#################################################################
# View model coefficients
#################################################################
resultsNames(dds)

#################################################################
# Log fold change shrinkage for visualization and ranking
#################################################################
library("apeglm")

#seems to run, unsure how long it takes
#resLFC <- lfcShrink(dds, coef="cell.type_syncytiotrophoblast_vs_composite", type="apeglm")
#resLFC

#################################################################
# MA-Plot
#################################################################

#unclear
plotMA(res)

#################################################################
# Extracting transformed values
#################################################################

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE) #runtime error/or takes long time
head(assay(vsd), 3)

#################################################################
# Effects of transformation on the variance
#################################################################

# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

#################################################################
# Data quality assessment by sample clustering and visualization
#################################################################
library("pheatmap")

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("cell.type","is.Paired", "sample.ID")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#################################################################
# Heatmap of the sample-to-sample distances
#################################################################

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$cell.type, vsd$is.Paired, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#################################################################
# PCA plot
#################################################################

plotPCA(vsd, intgroup=c("cell.type", "is.Paired"))


#is.Paired decoded by shape, easier to read
pcaData <- plotPCA(vsd, intgroup=c("cell.type", "is.Paired"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=cell.type, shape=is.Paired)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


plotPCA(vsd, intgroup=c("cell.type", "sample.ID"))

library("ggplot2")
pcaData <- plotPCA(vsd, intgroup=c("cell.type", "sample.ID"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=cell.type, shape=sample.ID)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#################################################################
# p-values and adjusted p-values
#################################################################

#p values in ascending order
resOrdered <- res[order(res$pvalue),]
summary(res)

#how many adjusted p-values were < .1
sum(res$padj < 0.1, na.rm=TRUE)

#set alpha to .05 and repeat
res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)