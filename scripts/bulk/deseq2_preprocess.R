
#set working directory ON CLUSTER
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


#OFF CLUSTER#

####################################################################################
# Saving the data as an R dataframe, both the meta and count data
####################################################################################
setwd("E:/Google Drive/Placenta_Cell_Types_Kyle/code/R code/local/placenta_sort/dex")


#get gex counts
counts.df <- read.table("placenta_sort_counts.txt", header = T, sep = '\t')

#use a featureCounts output file to get a vector of gene names
featureCounts <- read.table("E:/Google Drive/Placenta_Cell_Types_Kyle/output/run_2285_featureCounts/Sample_105812/Sample_105812+=_featureCounts.txt",
                            header = T,
                            stringsAsFactors = F,
                            sep = '\t')
gene_names = as.vector(featureCounts$Geneid, mode = "character")
gene_annotation = as.data.frame(featureCounts[1:6])

#add gene names as row names for the count matrix
counts.df <- as.data.frame(counts.df, row.names = gene_names)

#get sample meta data
meta.df <- read.table("meta_data.txt", header = T, sep = '\t', stringsAsFactors = F)

#convert meta data to dataframe
meta.df <- as.data.frame(meta.df)

#dropping data directory information
meta.df <- meta.df[,-5]

#saving the dataframes as r objects to preempt need to prepare the data again
saveRDS(counts.df, file = "sample_gex_counts_df.Rda")
saveRDS(meta.df, file = "sample_meta_data_df.Rda")
saveRDS(gene_annotation, file = "edgeR_gene_annotation.Rda")

test_load_counts <- readRDS("sample_gex_counts_df.Rda")
test_meta_data <- readRDS("sample_meta_data_df.Rda")

#code to rename extravillous trophoblast -> extravilloustrophoblast
levels(meta$cell.type)[levels(meta$cell.type) == "extravillous trophoblast"] <- "extravilloustrophoblast"


#counts are gex counts for each of the sample with the row names being the gene names
#meta stores the sample meta data
counts_uncoll <- readRDS("DESeq2_uncollapsed_count_data.Rda")
meta_uncoll <- readRDS("DESeq2_uncollapsed_sample_metadata.Rda")

#code to rename extravillous trophoblast -> extravilloustrophoblast
levels(meta_uncoll$cell.type)[levels(meta_uncoll$cell.type) == "extravillous trophoblast"] <- "extravilloustrophoblast"

counts <- readRDS(file = "DESeq2_collapsed_count_data.Rda")
meta <- readRDS(file = "DESeq2_collapsed_sample_metadata.Rda")

#code to rename extravillous trophoblast -> extravilloustrophoblast
levels(meta$cell.type)[levels(meta$cell.type) == "extravillous trophoblast"] <- "extravilloustrophoblast"

#re-saving cleaned data

saveRDS(meta_uncoll, file = "edgeR_uncollapsed_metadata.Rda")
saveRDS(meta, file = "edgeR_collapsed_metadata.Rda")
saveRDS(counts_uncoll, file = "edgeR_uncollapsed_counts.Rda")
saveRDS(counts, file = "edgeR_collapsed_counts.Rda")

####################################################################################
####################################################################################
# Reading in the clean data
####################################################################################
####################################################################################
setwd("E:/Google Drive/Placenta_Cell_Types_Kyle/code/R code/local/placenta_sort/dex/data/analysis_datasets")


#counts are gex counts for each of the sample with the row names being the gene names
#meta stores the sample meta data
counts_uncoll <- readRDS("edgeR_uncollapsed_counts.Rda")
meta_uncoll <- readRDS("edgeR_uncollapsed_metadata.Rda")

counts <- readRDS(file = "edgeR_collapsed_counts.Rda")
meta <- readRDS(file = "edgeR_collapsed_metadata.Rda")


gene_annotation <- readRDS("edgeR_gene_annotation.Rda")

#counts <- readRDS(file =   "./data/analysis_datasets/edgeR_collapsed_counts.Rda")
#meta <- readRDS(file =     "./data/analysis_datasets/edgeR_collapsed_metadata.Rda")
#gene_annotation <- readRDS("./data/analysis_datasets/edgeR_gene_annotation_with_symbols.Rda")

# Must make meta rownames identical to counts colnames
#rownames(meta) <- substr(rownames(meta), 8, 14)

# Creating the DESeq2 experiment object from count data and meta info
#dds <- DESeqDataSetFromMatrix(countData = counts,
#                              colData = meta,
#                              design = ~ 0 + sample.ID + cell.type)

# Adding comprehensive gene annotation to dds
#mcols(dds) <- DataFrame(mcols(dds), gene_annotation)

# Minimal pre-processing to reduce computational load - does NOT represent meaningful filtering
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

# Fitting the model
#dds <- DESeq(dds)

# Chop off the .xx from the ENSEMBL ID
#row.names(dds) <- substr(x = row.names(dds), start = 1, stop = 15)

#row.names(dds) <- mapIds(org.Hs.eg.db,
#                     keys=row.names(dds),
#                     column="SYMBOL",
#                     keytype="ENSEMBL",
#                     multiVals="first")


#BiomaRt returns 58,069 of 58,381 gene features from DESeq2 results. Of the 58,069 results, 56,521 return unique external gene names and 58,024 return unique ENSEMBL gene IDs.

# Setting up biomaRt
# Assigning Mart of interest, human ensembl here
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# Querying biomaRt using gene names to get back ENSEMBL IDs
DESeq2GeneIDs <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                       filters = c('ensembl_gene_id'),
                       values = rownames(dds),
                       mart = mart)
genes.tbl <- as_tibble(DESeq2GeneIDs)

#save/load
#saveRDS(DESeq2GeneIDs_tbl, "biomaRt_DESeq2_ensembl_ids.rda")
#DESeq2GeneIDs_tbl <- as_tibble(readRDS("biomaRt_DESeq2_ensembl_ids.rda"))

# Saved the fitted model object to avoid having to reestimate parameters
#saveRDS(dds, file = "./data/analysis_datasets/DESeq2_dds_noint_celltype_sampleid_nofiltering.rda")
#dds <- readRDS(file = "./data/analysis_datasets/DESeq2_dds_noint_celltype_sampleid_nofiltering.rda")

# Load the finished DESeq2 object - note this object has trimmed ensembl IDs
#dds <- readRDS("C:/Users/Kyle/Google Drive/Placenta_Cell_Types_Kyle/code/R code/local/placenta_sort/dex/data/analysis_datasets/DESeq2_dds_noint_celltype_sampleid_nofiltering.rda")
