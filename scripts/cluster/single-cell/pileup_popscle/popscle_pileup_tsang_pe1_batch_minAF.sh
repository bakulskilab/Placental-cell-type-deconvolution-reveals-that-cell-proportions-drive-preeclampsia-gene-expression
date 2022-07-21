#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=popscle_pileup_tsang_pe1_minAF
#SBATCH --mail-user=kyleac@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=160g
#SBATCH --time 168:00:00
#SBATCH --account=bakulski1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/scripts/pileup_popscle/%x-%j.log

# 12/6/20 - change to one node with 128 gigabytes
# 12/7/20 - use new vcf reference based on github.com/statgen/popscle/issues/25, see 1000_genomes_reference.sh in ~placenta_single_cell/reference
# 12/11/20 - updated to 1 week walltime

# make output dir
rm -r /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/tsang_pe1_minAF
mkdir /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/tsang_pe1_minAF

# unzip barcodes matrix if necessary
gunzip /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/tsang_pe1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

# run the popscle command through install on turbo
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle/bin/popscle dsc-pileup \
--sam /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/prep_popscle_pileup/tsang_pe1.minAF.bam \
--vcf /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
--out /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/tsang_pe1_minAF/pileupOUT \
--group-list /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/tsang_pe1/outs/filtered_feature_bc_matrix/barcodes.tsv

