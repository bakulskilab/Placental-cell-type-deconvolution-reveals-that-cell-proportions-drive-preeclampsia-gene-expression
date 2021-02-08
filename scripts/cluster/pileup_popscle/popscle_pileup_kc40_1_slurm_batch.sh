#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=popscle_pileup_kc40_1
#SBATCH --mail-user=kyleac@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=128000m
#SBATCH --time 12:00:00
#SBATCH --account=bakulski1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/scripts/pileup_popscle/%x-%j.log

# make output dir
rm -r /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/kc40_1
mkdir /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/kc40_1

# run the popscle command through install on turbo with CellRanger filtered cell barcodes
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle/bin/popscle dsc-pileup \
--sam /nfs/turbo/bakulski1/Datasets/Placenta_Single_Cell/ruddle.brcf.med.umich.edu/Run_2250/KC40_1_placenta/outs/possorted_genome_bam.bam \
--vcf /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/1000g_ref/1000g_ref_sorted_as_in_bam.vcf \
--out /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/kc40_1/dsc_pileup.kc40_1.pooled \
--group-list /nfs/turbo/bakulski1/Datasets/Placenta_Single_Cell/ruddle.brcf.med.umich.edu/Run_2250/KC40_1_placenta/outs/filtered_gene_bc_matrices/hg19/barcodes.tsv