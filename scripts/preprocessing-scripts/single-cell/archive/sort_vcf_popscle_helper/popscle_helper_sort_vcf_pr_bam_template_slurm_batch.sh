#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=popscle_helper_sortVCF_pr_bam_template_1000g_ref
#SBATCH --mail-user=kyleac@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=6GB
#SBATCH --time 12:00:00
#SBATCH --account=bakulski1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/1000g_ref/%x-%j.log

module load Bioinformatics
module load samtools
module load bcftools

dir="/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/1000g_ref/"

# Give execute permission to script
chmod +x $dir/sort_vcf_same_as_bam.sh

# Run popscle_helper script to sort 1000g_ref .vcf according to the same chromosome order in 10x output .bam
$dir/sort_vcf_same_as_bam.sh \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/possorted_genome_bam.bam \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/1000g_ref/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf \
> $dir/1000g_ref_sorted_as_in_pr_bam.vcf