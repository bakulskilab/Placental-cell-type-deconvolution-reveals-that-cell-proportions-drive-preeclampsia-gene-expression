#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=freemuxlet_pr478
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

rm -r /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478/results
mkdir /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478/results

# run the freemuxlet command through install on turbo
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle/bin/popscle freemuxlet \
--plp /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478/pileupOUT \
--out /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/pr_478/results/freemuxlet \
--nsample 2