G:\My Drive\Placenta_Cell_Types\RNA\placenta_cell_types_rna\scripts

Local analysis scripts

Sorted cell type bulk RNA-seq analysis:
1. placenta_sort_DESeq2_all_cell_types_prefilter
	input: featureCounts output

Single-cell analysis:
1. freemuxlet_analysis
	Summarizing freemuxlet results to be used for annotation of single-cell RNA-seq data; not necessary for downstream analysis
	input: HPC freemuxlet results
	
2. cellRanger_4.0_filtered_scRNAseq_qc
	QC of single-cell RNA-seq data
	input: HPC CellRanger results, (1)
	
3. seurat_annotation
	Annotation of single-cell RNA-seq data and exploratory clustering
	input: (2)
	
4. seurat_fetal_cell_type_assignment_biorep_mnn
	Fine-scale doublet removal, clustering, and annotation of fetal cell types. These results were used to annotate newly added Tsang et al. (2017) dataset.
	input: (3)
	
5. seurat_maternal_cell_type_assignment_biorep_mnn
	Fine-scale doublet removal, clustering, and annotation of maternal cell types. These results were used to annotate newly added Tsang et al. (2017) dataset.
	input: (3),(4)
	
6. seurat_dex
	Summary of cluster-specific (one vs. all) differential expression results and ontology
	input: (5)

7. seurat_analysis
	All other single-cell analyses including final QC results, cell type-specific gene expression, technical correlation, comparing proliferative vs. non-proliferative cytotrophoblasts, and various visualizations
	input: (5)

Comparing single-cell and sorted-cell expression results
1. bulk_single_cell_expression_correlation
	Validation of single-cell results using sorted bulk RNA-seq results
	input: (Sorted cell type bulk analysis) and (Single-cell analysis 5)

Deconvolution and Secondary Data Analysis of GSE75010
1. GSE75010_cell_type_proportion prep
	Download and prepare GSE75010 dataset and CIBERSORTx deconvolution results for analysis
	input: (CIBERSORTx GSE75010 deconvolution output)

2. GSE75010_cell_type_proportion_beta_regression_no_min_expr
	Cell type proportion PCA and beta regression analysis
	input: (1)

3. GSE75010_dex_analysis
	Gene expression PCA, limma differential expression, preparation of GSEA input, and mediation analysis
	input: (2)

4. GSE75010_gsea_res
	Analysis of differential expression GSEA results
	input: (GSEA 3 output)

5. GSE75010_mediation
	Mediation analysis with preeclampsia as the exposure, cell type composition as the mediator, and gene expression as the outcome
	input: (2), (3)

Deconvolution validation analyses
1. deconvolution_validation
	In silico testing of deconvolution reference and comparision of healthy whole tissue cell type proportions

2. facs_mixtures_deconvolution_analysis
	Various analyses of the Sorted cell type samples