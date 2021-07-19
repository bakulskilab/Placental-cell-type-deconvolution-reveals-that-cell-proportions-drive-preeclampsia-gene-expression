library(dplyr)
library(Seurat)
library(cowplot)
library(rgl)
library(EnvStats)
library(ggplot2)

# Working in the downloads folder temporarily
setwd("C:/Users/k_cam/Downloads/scRNA_hpc_qc")

#function to load in file, with sequential filtering steps
#first (droplets): no filtering
#second (filter_1000gene): 1000 gene minimum
#third (filter_3cell): min cell 3
filteringSteps <- function(file.name,obj.name)
{
  data.raw <- Read10X(data.dir = paste0("/nfs/turbo/bakulski1/Datasets/Placenta_Single_Cell/V2/ruddle.brcf.med.umich.edu/Run_2250/", file.name, "/outs/raw_gene_bc_matrices/hg19/"))
  #data.raw <- Read10X(data.dir=paste0("/nfs/turbo/bakulski1/Datasets/Pb_Hipp/ruddle.brcf.med.umich.edu/CellRanger_Output/",file.name,"_hippo/outs/raw_gene_bc_matrices/mm10/"))
  
  droplets <- CreateSeuratObject(counts = data.raw,min.features=0,min.cells=0,project = obj.name)
  
  filter_1000gene <- CreateSeuratObject(counts = data.raw,min.features=1000,min.cells=0,project = obj.name)
  
  filter_3cell <- CreateSeuratObject(counts = data.raw,min.features=1000,min.cells=3,project = obj.name)
  
  return(list(droplets=droplets,filter_1000gene=filter_1000gene,filter_3cell=filter_3cell))
}

#Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
# Warning when using the filter function
filter <- function(data.raw, obj.name) {
  
  droplets <- CreateSeuratObject(counts = data.raw,min.features=0,min.cells=0,project = obj.name)
  
  filter_1000gene <- CreateSeuratObject(counts = data.raw,min.features=1000,min.cells=0,project = obj.name)
  
  filter_3cell <- CreateSeuratObject(counts = data.raw,min.features=1000,min.cells=10,project = obj.name)  

  return(list(droplets=droplets,filter_1000gene=filter_1000gene,filter_3cell=filter_3cell))
}

# Functions that uses the add 1 minus 1 workaround for zeros in the geometric mean or geometric standard deviation function from EnvStats
# Input is a non-negative vector
# Checks that zeros are in the column and uses the default mode if no zeros are detected
# "Geometric mean extension for data sets with zeros", arXiv:1806.06403 helpful for a short summary on GeoMean and zeros
# http://alexkritchevsky.com/2018/06/15/geometric-mean.html Another Geometric stats resource
geo_Mean_zeros <- function(data) {
  if(0 %in% data) {
    data <- data+1
    result <- geoMean(data)
    result <- result - 1
    return(result)
  }
  else {
    return(geoMean(data))
  }
}

# TODO: appropriate workaround for the geometric standard deviation below, not sure if the result - 1 workaround is appropriate for geoSD
geo_SD_zeros <- function(data) {
  if(0 %in% data) {
    data <- data+1
    result <- geoSD(data)
    return(result)
  }
  else {
    return(geoSD(data))
  }
}

#TODO: pctMT doesn't work for droplets
QCValues_zeros <- function(SeuratList){
  QC.values <- data.frame(matrix(0,nrow=3,ncol=8))
  rownames(QC.values) <- c('NoFilter','Filter1000g','Filter1000g_3cell')
  colnames(QC.values) <- c('Droplets','Mito_Gt_15pct','gm_pctMT','gsd_pctMT','gm_nCount','gsd_nCount','gm_nFeature','gsd_nFeature')
  
  pct.mito.drop <- PercentageFeatureSet(object = SeuratList$droplets, pattern = "^MT-")
  SeuratList$droplets <- AddMetaData(object = SeuratList$droplets, metadata = pct.mito.drop, col.name = "percent.mito")
  
  pct.mito.1000g <- PercentageFeatureSet(object = SeuratList$filter_1000gene, pattern = "^MT-")
  SeuratList$filter_1000gene <- AddMetaData(object = SeuratList$filter_1000gene, metadata = pct.mito.1000g, col.name = "percent.mito")
  
  pct.mito.3c <- PercentageFeatureSet(object = SeuratList$filter_3cell, pattern = "^MT-")
  SeuratList$filter_3cell <- AddMetaData(object = SeuratList$filter_3cell, metadata = pct.mito.3c, col.name = "percent.mito")
  
  QC.values$Droplets <- c(table(SeuratList$droplets@meta.data$nCount_RNA>0)['TRUE'],
                          table(SeuratList$filter_1000gene@meta.data$nCount_RNA>0)['TRUE'],
                          table(SeuratList$filter_3cell@meta.data$nCount_RNA>0)['TRUE'])
  QC.values$Mito_Gt_15pct <- c(length(pct.mito.drop[pct.mito.drop>15 & !is.na(pct.mito.drop)]),
                               length(pct.mito.1000g[pct.mito.1000g>15 & !is.na(pct.mito.1000g)]),
                               length(pct.mito.3c[pct.mito.3c>15 & !is.na(pct.mito.3c)]))
  QC.values$gm_pctMT <- c(geo_Mean_zeros(pct.mito.drop[,1]),
                          geo_Mean_zeros(pct.mito.1000g[,1]),
                          geo_Mean_zeros(pct.mito.3c[,1]))
  QC.values$gsd_pctMT <- c(geo_SD_zeros(pct.mito.drop[,1]),
                           geo_SD_zeros(pct.mito.1000g[,1]),
                           geo_SD_zeros(pct.mito.3c[,1]))
  QC.values$gm_nCount <- c(geo_Mean_zeros(SeuratList$droplets@meta.data$nCount_RNA),
                           geo_Mean_zeros(SeuratList$filter_1000gene@meta.data$nCount_RNA),
                           geo_Mean_zeros(SeuratList$filter_3cell@meta.data$nCount_RNA))
  QC.values$gsd_nCount <- c(geo_SD_zeros(SeuratList$droplets@meta.data$nCount_RNA),
                            geo_SD_zeros(SeuratList$filter_1000gene@meta.data$nCount_RNA),
                            geo_SD_zeros(SeuratList$filter_3cell@meta.data$nCount_RNA))
  QC.values$gm_nFeature <- c(geo_Mean_zeros(SeuratList$droplets@meta.data$nFeature_RNA),
                             geo_Mean_zeros(SeuratList$filter_1000gene@meta.data$nFeature_RNA),
                             geo_Mean_zeros(SeuratList$filter_3cell@meta.data$nFeature_RNA))
  QC.values$gsd_nFeature <- c(geo_SD_zeros(SeuratList$droplets@meta.data$nFeature_RNA),
                              geo_SD_zeros(SeuratList$filter_1000gene@meta.data$nFeature_RNA),
                              geo_SD_zeros(SeuratList$filter_3cell@meta.data$nFeature_RNA))
  
  SeuratList$QC <- QC.values
  return(SeuratList)
}


# Looks somewhat lognormally distributed
ggplot(data = pct.mito.1000g.log) + aes(nCount_RNA) + geom_density()

# TODO how to do this qqplot as another way to check lognormality
ggplot(pct.mito.1000g.log) + geom_qq_line(aes(sample = nCount_RNA))

# Read in the raw data and use filter to filter with default thresholds used for KC40 & KC42 (1000 genes and 3 cells)
tnl.478.raw <- readRDS("tnl_478_raw.rda")
tnl.481.raw <- readRDS("tnl_481_raw.rda")
tnl.484.raw <- readRDS("tnl_484_raw.rda")

tnl.478 <- filter(tnl.478.raw, "tnl.478")
tnl.481 <- filter(tnl.481.raw, "tnl.481")
tnl.484 <- filter(tnl.484.raw, "tnl.484")

tnl.478 <- QCValues_zeros(tnl.478)
tnl.481 <- QCValues_zeros(tnl.481)
tnl.484 <- QCValues_zeros(tnl.484)

tnl.list <- list(tnl.478 = tnl.478,
                  tnl.481 = tnl.481,
                  tnl.484 = tnl.484)

saveRDS(tnl.list, "tnl_list.rda")

kc40.1.raw <- readRDS(file = "kc40_1_raw.rda")
kc40.2.raw <- readRDS(file = "kc40_2_raw.rda")
kc42.1.raw <- readRDS(file = "kc42_1_raw.rda")
kc42.2.raw <- readRDS(file = "kc42_2_raw.rda")

kc40.1 <- filter(kc40.1.raw, "kc.40.1")
kc40.2 <- filter(kc40.2.raw, "kc.40.2")
kc42.1 <- filter(kc42.1.raw, "kc.42.1")
kc42.2 <- filter(kc42.2.raw, "kc.42.2")

kc40.1 <- QCValues_zeros(kc40.1)
kc40.2 <- QCValues_zeros(kc40.2)
kc42.1 <- QCValues_zeros(kc42.1)
kc42.2 <- QCValues_zeros(kc42.2)

kc.list <- list(kc40.1 = kc40.1,
                kc40.2 = kc40.2,
                kc42.1 = kc42.1,
                kc42.2 = kc42.2)

saveRDS(kc.list, "kc_list.rda")


