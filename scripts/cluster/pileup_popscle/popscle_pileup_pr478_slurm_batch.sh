#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=popscle_pileup_pr478
#SBATCH --mail-user=kyleac@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=128000m
#SBATCH --time 12:00:00
#SBATCH --account=bakulski1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/%x-%j.log

# 12/6/20 - change to one node with 128 gigabytes

# make output dir
rm -r /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478
mkdir /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478

# unzip barcodes matrix
#gunzip /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

# run the popscle command through install on turbo
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle/bin/popscle dsc-pileup \
--sam /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/possorted_genome_bam.bam \
--vcf /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/1000g_ref/1000g_ref_sorted_as_in_bam.vcf \
--out /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478/pr_478_pooled_dsc_pileup \
--group-list /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/filtered_feature_bc_matrix/barcodes.tsv

# re-zip barcodes matrix, untested
#gzip /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/filtered_feature_bc_matrix/barcodes.tsv