#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=batch_pique-regi_tnl_481
#SBATCH --mail-user=kyleac@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=6GB
#SBATCH --time 48:00:00
#SBATCH --account=bakulski1
#SBATCH --partition=standard
#SBATCH --output=/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/%x-%j.log

/scratch/bakulski_root/bakulski1/kyleac/cellranger/app/cellranger-4.0.0/bin/cellranger \
count \
--id=batch_pique_regi_tnl_481 \
--fastqs=/scratch/bakulski_root/bakulski1/kyleac/dbgap/fastq \
--sample=SRR10166481 \
--transcriptome=/scratch/bakulski_root/bakulski1/kyleac/cellranger/human_reference/refdata-gex-GRCh38-2020-A