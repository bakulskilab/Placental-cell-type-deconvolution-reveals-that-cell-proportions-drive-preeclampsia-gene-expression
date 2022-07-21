#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=bam2fastq_v1.2.1_tsang_PN1
#SBATCH --mail-user=kyleac@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=6GB
#SBATCH --time 48:00:00
#SBATCH --account=bakulski1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/%x-%j.log

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/lib/bin/bamtofastq --cr11 \
/nfs/turbo/bakulski1/People/kyleac/tsang_2017/tsang_scRNA_bams/PN1/PN1.bam \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PN1
