#running both Placenta Sort experiments on the same multiqc
#Sample105* were paired end whereas Sample103* were single end
#runs locally

#load python
module load python-anaconda3

#execute python script that installs multiqc
#import_multiqc.py
python import_multiqc.py

multiqc \
/nfs/turbo/bakulski1/People/kyleac/Placenta_Sort_RNA/Run_2245_output/Sample_103* \
/nfs/turbo/bakulski1/People/kyleac/Placenta_Sort_RNA/run_2285_output/fastQC/Sample_105*