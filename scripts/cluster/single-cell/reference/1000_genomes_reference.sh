# https://arc-ts.umich.edu/wp-content/uploads/sites/4/2020/05/Great-Lakes-Cheat-Sheet.pdf GreatLakes cheat sheet
# https://github.com/SchlossLab/Great_Lakes_SLURM slurm guide

# OUTPUT:
# 1000g.minAF-0.1.reheader.vcf.gz used for freemuxlet reference file
# Each sequencing bam file with MAF < 10% filtered out (sampleID.minAF.bam) used for freemuxlet input

# December 11th, 2020
# Potential MUCH faster alternative from Aerts lab member: https://github.com/aertslab/popscle_helper_tools/issues/4
# Download just the SNV data from 1000g and filter out the rare variants, reduces file sizes by orders of magnitude and
# filters out >85% of total sites (rare variant sites (MAF < 10%)) since rare variants <5% MAF vastly outnumber >% MAF
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz
# No index, add with bcftools
bcftools index ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz
# Filter out rare variants, reference doesn't have Y or MT
bcftools view -r \
    chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX \
    -i 'INFO/AF[0] < 0.90 && INFO/AF[0] > 0.10' \
    ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz \
    -O z -o 1000Genomes.wgs.GRCH38.sites.minAF-0.1.freemuxlet.vcf.gz
	
# Compare .bam vs. .vcf
# Do not differ by CHROM naming convention, i.e. chr1 vs. chr1
# Contigs ordered differently - lexico vs. numeric

# https://www.biostars.org/p/133487/ Pierre Lindenbaum's solution to sort by chr lexicographically and then by position
grep '^#' 1000Genomes.wgs.GRCH38.sites.minAF-0.1.freemuxlet.vcf > 1000g.minAF-0.1.out.vcf \
 && grep -v '^#' 1000Genomes.wgs.GRCH38.sites.minAF-0.1.freemuxlet.vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> 1000g.minAF-0.1.out.vcf

# Bzip for bcftools input and index
bgzip -c 1000g.minAF-0.1.out.vcf > 1000g.minAF-0.1.out.vcf.gz # compress
bcftools index 1000g.minAF-0.1.out.vcf.gz					  # index for bcftools operations

# Download header.txt and subset header (with Notepad) to well-defined regions outlined (eliminate all contigs besides 1-22, X, Y as above and then reorder lexicographically
# Upload subset and lexicographically sorted header.txt back to the cluster
# Reheader old file with manually curated new header
bcftools view 1000g.minAF-0.1.out.vcf.gz -h > minAF.header.txt

# After manual notepad editing of header, reheader reference file
bcftools reheader -h minAF.header.txt -o 1000g.minAF-0.1.reheader.vcf.gz 1000g.minAF-0.1.out.vcf.gz

# Add relevant dependences to PATH
PATH=$PATH:/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/dependencies

# Filter .bams for all samples because so effective at speeding up pileup
## PR Samples
./filter_bam_file_for_popscle_dsc_pileup.sh \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/possorted_genome_bam.bam \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/filtered_feature_bc_matrix/barcodes.tsv \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
./pr478.minAF.bam

./filter_bam_file_for_popscle_dsc_pileup.sh \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_481/outs/possorted_genome_bam.bam \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_481/outs/filtered_feature_bc_matrix/barcodes.tsv \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
./pr481.minAF.bam

./filter_bam_file_for_popscle_dsc_pileup.sh \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_484/outs/possorted_genome_bam.bam \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_484/outs/filtered_feature_bc_matrix/barcodes.tsv \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
./pr484.minAF.bam

# KC samples
### RUN this after CELLRANGER RUN
./filter_bam_file_for_popscle_dsc_pileup.sh \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_40_1/outs/possorted_genome_bam.bam \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_40_1/outs/filtered_feature_bc_matrix/barcodes.tsv \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
./kc40.1.minAF.bam

./filter_bam_file_for_popscle_dsc_pileup.sh \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_40_2/outs/possorted_genome_bam.bam \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_40_2/outs/filtered_feature_bc_matrix/barcodes.tsv \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
./kc40.2.minAF.bam

./filter_bam_file_for_popscle_dsc_pileup.sh \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_42_1/outs/possorted_genome_bam.bam \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_42_1/outs/filtered_feature_bc_matrix/barcodes.tsv \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
./kc42.1.minAF.bam

./filter_bam_file_for_popscle_dsc_pileup.sh \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_42_2/outs/possorted_genome_bam.bam \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_42_2/outs/filtered_feature_bc_matrix/barcodes.tsv \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
./kc42.2.minAF.bam

# Tsang samples
# 1/31/22 
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/scripts/popscle_helper_tools-master/filter_bam_file_for_popscle_dsc_pileup.sh \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/tsang_pe1/outs/possorted_genome_bam.bam \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/tsang_pe1/outs/filtered_feature_bc_matrix/barcodes.tsv \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/1000g.minAF-0.1.reheader.vcf \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/prep_popscle_pileup/tsang_pe1.minAF.bam
