#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=freemuxlet_batch_minAF
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

for path in /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/*minAF; do

	SAMPLEID="$(basename "${path}")"
	
	rm -r /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$SAMPLEID/results
	mkdir /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$SAMPLEID/results

	/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle/bin/popscle freemuxlet \
	--plp /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$SAMPLEID/pileupOUT \
	--out /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$SAMPLEID/results/$SAMPLEID \
	--nsample 2

	FILENAME="${SAMPLEID}.clust1.samples.gz"	
	OUTPUTNAME="${SAMPLEID}.clusters"
	
	zcat /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$SAMPLEID/results/$FILENAME \
	> /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$SAMPLEID/results/$OUTPUTNAME
	
done