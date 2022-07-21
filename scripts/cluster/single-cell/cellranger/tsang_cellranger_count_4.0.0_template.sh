#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=cellranger_count_tsang2017
#SBATCH --mail-user=kyleac@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=6GB
#SBATCH --time 168:00:00
#SBATCH --account=bakulski1
#SBATCH --partition=standard
#SBATCH --output=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/%x-%j.log

cd /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/bin/cellranger \
count \
--id=tsang_pn1 \
--fastqs=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PN1/gemgroup001 \
--transcriptome=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/old_cellranger/human_reference/refdata-gex-GRCh38-2020-A \
--chemistry=SC3Pv1

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/bin/cellranger \
count \
--id=tsang_pn2 \
--fastqs=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PN2/gemgroup001 \
--transcriptome=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/old_cellranger/human_reference/refdata-gex-GRCh38-2020-A \
--chemistry=SC3Pv1

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/bin/cellranger \
count \
--id=tsang_pn3c \
--fastqs=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PN3C/gemgroup001 \
--transcriptome=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/old_cellranger/human_reference/refdata-gex-GRCh38-2020-A \
--chemistry=SC3Pv1

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/bin/cellranger \
count \
--id=tsang_pn3p \
--fastqs=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PN3P/gemgroup001 \
--transcriptome=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/old_cellranger/human_reference/refdata-gex-GRCh38-2020-A \
--chemistry=SC3Pv1

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/bin/cellranger \
count \
--id=tsang_pn4c \
--fastqs=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PN4C/gemgroup001 \
--transcriptome=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/old_cellranger/human_reference/refdata-gex-GRCh38-2020-A \
--chemistry=SC3Pv1

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/bin/cellranger \
count \
--id=tsang_pn4p \
--fastqs=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PN4P/gemgroup001 \
--transcriptome=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/old_cellranger/human_reference/refdata-gex-GRCh38-2020-A \
--chemistry=SC3Pv1

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/bin/cellranger \
count \
--id=tsang_pe1 \
--fastqs=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PE1/gemgroup001 \
--transcriptome=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/old_cellranger/human_reference/refdata-gex-GRCh38-2020-A \
--chemistry=SC3Pv1

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/bin/cellranger \
count \
--id=tsang_pe2 \
--fastqs=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PE2/gemgroup001 \
--transcriptome=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/old_cellranger/human_reference/refdata-gex-GRCh38-2020-A \
--chemistry=SC3Pv1

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/bin/cellranger \
count \
--id=tsang_pe3 \
--fastqs=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PE3/gemgroup001 \
--transcriptome=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/old_cellranger/human_reference/refdata-gex-GRCh38-2020-A \
--chemistry=SC3Pv1

/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/cellranger-4.0.0/bin/cellranger \
count \
--id=tsang_pe4 \
--fastqs=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/bamtofastq_output/PE4/gemgroup001 \
--transcriptome=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/old_cellranger/human_reference/refdata-gex-GRCh38-2020-A \
--chemistry=SC3Pv1