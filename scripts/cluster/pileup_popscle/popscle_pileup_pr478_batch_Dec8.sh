#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=popscle_pileup_pr478
#SBATCH --mail-user=kyleac@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=160g
#SBATCH --time 12:00:00
#SBATCH --account=bakulski1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/scripts/pileup_popscle/%x-%j.log

# 12/6/20 - change to one node with 128 gigabytes
# 12/7/20 - use new vcf reference based on github.com/statgen/popscle/issues/25, see 1000_genomes_reference.sh in ~placenta_single_cell/reference

# make output dir
rm -r /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478
mkdir /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478

# unzip barcodes matrix if necessary
gunzip /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

# run the popscle command through install on turbo
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle/bin/popscle dsc-pileup \
--sam /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/possorted_genome_bam.bam \
--vcf /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/with_chr.vcf \
--out /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478/pileupOUT \
--group-list /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/filtered_feature_bc_matrix/barcodes.tsv