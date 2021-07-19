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

# Background
# https://github.com/statgen/popscle/issues/9 explains that freemuxlet reference file is to get AFs

# December 8th, trying to sort vcf lexicographically and replace add chr prefix to 1000g reference
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