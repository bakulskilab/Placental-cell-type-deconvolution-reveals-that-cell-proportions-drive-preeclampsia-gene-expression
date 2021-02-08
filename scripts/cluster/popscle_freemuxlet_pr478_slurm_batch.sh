#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=popscle_freemuxlet_pr478
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

# run the popscle command through install on turbo
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle/bin/popscle freemuxlet \
--plp /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478/pr_478_pooled_dsc_pileup \
--nsample 2 \
--out /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478/pr_478_freemuxlet_pooled
