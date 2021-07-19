## Align to genome and infer libtype
# https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype/

# Load Bionformatics tools
module load Bioinformatics
module load salmon

# Salmon single-end command
/nfs/turbo/bakulski1/People/kyleac/Placenta_Sort_RNA/salmon-1.5.1_linux_x86_64/bin/salmon quant \
-i /nfs/turbo/bakulski1/People/kyleac/Placenta_Sort_RNA/salmon_ref/alias/hg38/salmon_partial_sa_index/default/ \
-l A \
-r /nfs/turbo/bakulski1/Datasets/Placenta_Sort_RNA/ruddle.brcf.med.umich.edu/Run_2245/colacino/Sample_103893/103893_ATTCAGAA-CCTATCCT_S1_L007_R1_001.fastq.gz \
-o /nfs/turbo/bakulski1/People/kyleac/Placenta_Sort_RNA/salmon_out/sample_103893 \
--validateMappings

# Confirm detected as SR (stranded - reverse)


# Salmon paired-end command
/nfs/turbo/bakulski1/People/kyleac/Placenta_Sort_RNA/salmon-1.5.1_linux_x86_64/bin/salmon quant \
-i /nfs/turbo/bakulski1/People/kyleac/Placenta_Sort_RNA/salmon_ref/alias/hg38/salmon_partial_sa_index/default/ \
-l A \
-1 /nfs/turbo/bakulski1/Datasets/Placenta_Sort_RNA/ruddle.brcf.med.umich.edu/Run_2285/bakulski/Sample_105812/105812_ATTACTCG-AGGCTATA_S1_L004_R1_001.fastq.gz \
-2 /nfs/turbo/bakulski1/Datasets/Placenta_Sort_RNA/ruddle.brcf.med.umich.edu/Run_2285/bakulski/Sample_105812/105812_ATTACTCG-AGGCTATA_S1_L004_R2_001.fastq.gz \
-o /nfs/turbo/bakulski1/People/kyleac/Placenta_Sort_RNA/salmon_out/sample_103893 \
--validateMappings

# ISR most likely libarary type for paired-end sample