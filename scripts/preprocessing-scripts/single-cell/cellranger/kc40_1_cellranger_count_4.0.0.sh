#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=batch_kc_40_1
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
--id=batch_kc_40_1 \
--fastqs=/nfs/turbo/bakulski1/Datasets/Placenta_Single_Cell/V2/ruddle.brcf.med.umich.edu/Run_2250/outs/fastq_path/HTCVVBBXX/KC40_1_placenta \
--sample=KC40_1_placenta \
--transcriptome=/scratch/bakulski_root/bakulski1/kyleac/cellranger/human_reference/refdata-gex-GRCh38-2020-A