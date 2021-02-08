# https://arc-ts.umich.edu/wp-content/uploads/sites/4/2020/05/Great-Lakes-Cheat-Sheet.pdf GreatLakes cheat sheet
# https://github.com/SchlossLab/Great_Lakes_SLURM slurm guide

# December 7th, 2020
# https://github.com/statgen/popscle/issues/25 Following issues/25 to get reference VCF for freemuxlet

# Download 1000genomes 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
# The file is 1.8 GB

# SORT, THEN REHEADER IF NEEDED, NOT THE OTHER WAY AROUND
# Add samtools to PATH
export PATH=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/dependencies/bin:$PATH
# Add bcftools to PATH
export PATH=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/dependencies/bcftools:$PATH
# Add bedtools to PATH
export PATH=/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/dependencies/bedtools2/bin:$PATH

################### PIPELINE BELOW USED  #############################
# Background
# https://github.com/statgen/popscle/issues/9 explains that freemuxlet reference file is get AFs

# December 8th, trying to sort vcf lexico graphically and replace add chr prefix to 1000g reference
# 1000 genome reference is numerically sorted with no chr prefix

# https://github.com/statgen/popscle/issues/35 example code to pull SQ (contig equivalent in .bam) 
# First, compare 10 .bam output and 1000g_ref.vcf
samtols view -h ./possorted_genome_bam.bam | less
zcat ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz | less

# Compare .bam vs. .vcf
# Differ by CHROM - chr1 vs. 1
# Differ by MT - M vs. MT
# Contigs ordered differently - lexico vs. numeric

# https://www.biostars.org/p/133487/ Pierre Lindenbaum's solution to sort by chr lexicographically and then by position
grep '^#' ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf > out.vcf \
 && grep -v '^#' ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> out.vcf
# https://github.com/statgen/demuxlet/issues/36 reheader after lexicographically sorting
# bgzip sorted output for bcftools input
bgzip -c out.vcf > out.vcf.gz
# Below commands produces no index found error
# https://www.biostars.org/p/336800/ Index bcftools file for subequent subset command
bcftools index out.vcf.gz
# Subset to well defined regions
bcftools view out.vcf.gz --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y > subset.vcf.gz
# Output header
bcftools view subset.vcf.gz -h > header.txt
# Download header.txt and subset header to well-defined regions outlined (eliminate all contigs besides 1-22, X, Y as above and then reorder lexicographically
# Upload subset and lexicographically sorted header.txt back to the cluster
# Reheader old file with manually curated new header
bcftools reheader -h header.txt -o subset_reheader.vcf subset.vcf.gz
# Created uncompressed copy for awk
#bgzip -dc > subset.vcf
# https://www.biostars.org/p/98582/ Code written by John under Vivek's reply to add chr prefix from chromosome notation in vcf file contents and header
# MODIFIED BY ME: add '[0-9]+' in matching string to not add chr to GL contigs, though those should be gone now that I deleted those earlier
awk '{ 
        if($0 !~ /^#/) 
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)([0-9]+.*)/,m))
            print m[1]"chr"m[2];
        else print $0 
      }' subset_reheader.vcf > with_chr.vcf
# 1000g_ref.vcf should be ready for pileup input
# Appears to work but very slow

# Try filtering .bam to speed up several hundred times with popscle_helper
$ ./filter_bam_file_for_popscle_dsc_pileup.sh
Usage:   filter_bam_file_for_popscle_dsc_pileup input_bam_filename barcodes_tsv_filename vcf_filename output_bam_filename

./filter_bam_file_for_popscle_dsc_pileup.sh \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/possorted_genome_bam.bam \
/scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/filtered_feature_bc_matrix/barcodes.tsv \
/nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/reference/with_chr.vcf \
./out.bam

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
# Do not differ by CHROM - chr1 vs. chr1
# Contigs ordered differently - lexico vs. numeric

# https://www.biostars.org/p/133487/ Pierre Lindenbaum's solution to sort by chr lexicographically and then by position
grep '^#' 1000Genomes.wgs.GRCH38.sites.minAF-0.1.freemuxlet.vcf > 1000g.minAF-0.1.out.vcf \
 && grep -v '^#' 1000Genomes.wgs.GRCH38.sites.minAF-0.1.freemuxlet.vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> 1000g.minAF-0.1.out.vcf

# Bzip for bcftools input and index
bgzip -c 1000g.minAF-0.1.out.vcf > 1000g.minAF-0.1.out.vcf.gz
bcftools index 1000g.minAF-0.1.out.vcf.gz

# Download header.txt and subset header to well-defined regions outlined (eliminate all contigs besides 1-22, X, Y as above and then reorder lexicographically
# Upload subset and lexicographically sorted header.txt back to the cluster
# Reheader old file with manually curated new header
bcftools view 1000g.minAF-0.1.out.vcf.gz -h > minAF.header.txt
bcftools reheader -h minAF.header.txt -o 1000g.minAF-0.1.reheader.vcf.gz 1000g.minAF-0.1.out.vcf.gz

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
