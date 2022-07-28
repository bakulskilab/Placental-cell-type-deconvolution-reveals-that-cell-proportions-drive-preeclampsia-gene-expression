/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell

/cellranger/ contains the cellranger software, associated scripts, and output
/data/ contains key files required for scripts run on the cluster
/dependencies/ contains software required for cluster analyses, including freemuxlet and tools to manipulate genetic data files such as .vcf or .bam
/fastqc_output/ contains fastQC output to verify sequencing read lengths
/popscle/ contains the popscle software suite
/popscle_out/ contains output of popscle software
/reference/ contains reference genome information for use with popscle and raw RNA-seq data processed for popscle input
/scripts/ contains the scripts required to run freemuxlet and popscle tools
/tutorial_dataset/ contains tutorial dataset for popscle tools

Analyses conducted and order of operations:

Freemuxlet demultiplexing:
1. 1000_genomes_reference.sh
	Code to format 1000g reference and raw single-cell RNA-seq data for input into popscle-pileup
2. popscle_pileup
	Scripts found in ~/scripts/pileup_popscle/
3. popscle-freemuxlet
	Script found in ~/scripts/freemuxlet
4. Results exported for local analysis in "freemuxlet_analysis.Rmd"

Cellranger mapping and feature counting
1. Run cellranger count on raw reads
	Scripts and results found in ~/cellranger/run_cellranger_count/
2. Download 10x filtered matrices to local machine for downstream analysis in "cellRanger_4.0_filtered_scRNAseq_qc.Rmd"
	Filtered matrices can be found in ~cellranger/filtered_10x_output/
	


