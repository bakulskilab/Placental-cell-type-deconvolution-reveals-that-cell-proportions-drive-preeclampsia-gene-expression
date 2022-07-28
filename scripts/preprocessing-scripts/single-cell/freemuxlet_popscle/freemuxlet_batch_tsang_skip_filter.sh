#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=freemuxlet_tsang_batch_minAF_skip_filter
#SBATCH --mail-user=kyleac@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=64g
#SBATCH --time 168:00:00
#SBATCH --account=bakulski1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/scripts/freemuxlet/%x-%j.log

# Add relevant dependences to PATH to run filter_bam_file_for_popscle_dsc_pileup.sh 
PATH=$PATH:/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/dependencies/bin
PATH=$PATH:/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/dependencies/bedtools2/bin

# Get a vector of pathnames for which cellranger output files you want to run freemuxlet pipeline on:
# 1. Filter to common variants and cell barcodes
# 2. Run popscle-pileup
# 3. Run freemuxlet
for path in /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/tsang*; do

	# Create a sample ID to be used in all filenames
	SAMPLEID="$(basename "${path}")"
	
	### 1. Filter cellranger output to common variants and cell barcodes using filter_bam_file_for_popscle_dsc_pileup.sh ###
	# Make filtered minAF output dir
	MINAF_FOLDER_NAME="${SAMPLEID}_minAF"
	#rm -r /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME
	#mkdir /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME
	
	# Unzip barcodes matrix
	gunzip /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/$SAMPLEID/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
	
	# Run filter script
	MINAF_FILE_NAME="${SAMPLEID}_minAF.bam"
	#/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/scripts/popscle_helper_tools-master/filter_bam_file_for_popscle_dsc_pileup.sh \
	#/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/$SAMPLEID/outs/possorted_genome_bam.bam \
	#/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/$SAMPLEID/outs/filtered_feature_bc_matrix/barcodes.tsv \
	#/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
	#/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME/$MINAF_FILE_NAME
	
	### 2. Run popscle-pileup on minAF-filtered .bams ###

	# run the popscle command through install on turbo
	/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle/bin/popscle dsc-pileup \
	--sam /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME/$MINAF_FILE_NAME \
	--vcf /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
	--out /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME \
	--group-list /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/$SAMPLEID/outs/filtered_feature_bc_matrix/barcodes.tsv
	
	# Rezip barcodes matrix
	gzip /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/$SAMPLEID/outs/filtered_feature_bc_matrix/barcodes.tsv
	
	### 3. Run freemuxlet on popscle-pileup output ###
	rm -r /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME/results
	mkdir /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME/results

	/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle/bin/popscle freemuxlet \
	--plp /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME \
	--out /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME/results \
	--nsample 2

	FILENAME="${SAMPLEID}.clust1.samples.gz"	
	OUTPUTNAME="${SAMPLEID}.clusters"
	
	zcat /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$SAMPLEID/results/$FILENAME \
	> /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$SAMPLEID/results/$OUTPUTNAME
	
done
