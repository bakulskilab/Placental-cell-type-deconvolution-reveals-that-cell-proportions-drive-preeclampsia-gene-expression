source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("apeglm")
biocLite("vsn")
biocLite("pcaExplorer")

install.packages("pheatmap")

library("png")
library("QoRTs")
library("JunctionSeq")
library("DESeq2")
library("DEXSeq")
library("edgeR")
library("data.table")
library("qdapTools")

library("pheatmap")

#set working directory
setwd("/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/analysis_r/DESeq2/")

#set sample directory names
run2245_dir_prefix = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/run2245_output/featureCounts"
run2285_dir_prefix = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/run2285_output/featureCounts"

#testing that geneid list in the same between samples
test <- read.table(file = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/run2245_output/featureCounts/Sample_103902/Sample_103902+=_007_featureCounts.txt", 
                   sep = '\t',
                   header = TRUE,
                   stringsAsFactors = FALSE)
test2 <- read.table(file = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/run2285_output/featureCounts/Sample_105815/Sample_105815+=_featureCounts.txt", 
                    sep = '\t',
                    header = TRUE,
                    stringsAsFactors = FALSE)

#geneid tests TRUE for between 2245/2285 and different samples within 2245
identical(test$Geneid, test2$Geneid)


#############################################
# Building featureCounts expression matrix
#############################################

#read in placenta sort decoder that's necessary to run this block of code
placenta_sort_decoder_df <- read.table(file = "/scratch/bakulski_fluxod/kyleac/placenta_sort_RNA/analysis_r/placenta_sort_decoder.txt",
                                       header = T, sep = '\t', stringsAsFactors = F)


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

rownames(counts) <- counts$Geneid
counts <- counts[,-1]

saveRDS(counts, file = "placenta_sort_counts.Rda")

###################################################
# Entering the data into DESeq2 and Vignette
###################################################
cts <- readRDS(file = "placenta_sort_counts.Rda")

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

#lane.ID and/or is.Paired were linear combinations of each other and so couldn't run in the model
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ sample.ID + is.Paired + cell.type)

#export uncollapsed data for portability
#sample data:
export_meta_data <- dds@colData
export <- as.data.frame(export_meta_data)
saveRDS(export, file = "DESeq2_uncollapsed_sample_metadata.Rda")

#counts
cts_export <- as.data.frame(dds@assays$data$counts)
identical(rownames(cts_export), test$Geneid)
saveRDS(cts_export, file = "DESeq2_uncollapsed_count_data.Rda")

#################################################################
# Collapse Technical Replicates
#################################################################

dds <- collapseReplicates(dds, dds$run.sample.ID, renameCols = TRUE)

#exporting collapsed data for portability:
#sample data
export_meta_data <- dds@colData
export <- as.data.frame(export_meta_data)
export <- export[, -c(1,2, 5)]
saveRDS(export, file = "DESeq2_collapsed_sample_metadata.Rda")

#counts
cts_export <- as.data.frame(dds@assays$data$counts)
colnames(cts_export)[1:16] <- substr(x = colnames(cts_export)[1:16], start = 1, stop = 6)
identical(rownames(cts_export), test$Geneid)
saveRDS(cts_export, file = "DESeq2_collapsed_count_data.Rda")
