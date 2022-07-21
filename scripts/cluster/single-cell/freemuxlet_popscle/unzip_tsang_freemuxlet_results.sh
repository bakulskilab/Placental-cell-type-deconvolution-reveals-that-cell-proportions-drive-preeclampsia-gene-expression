for path in /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/data/cellranger_count_output/tsang*; do

	# Create a sample ID to be used in all filenames
	SAMPLEID="$(basename "${path}")"
	MINAF_FOLDER_NAME="${SAMPLEID}_minAF"
	FILENAME="${SAMPLEID}.clust1.samples.gz"	
	OUTPUTNAME="${SAMPLEID}.clusters"
	
	zcat /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME/$FILENAME \
	> /nfs/turbo/bakulski1/People/kyleac/placenta_single_cell/popscle_out/$MINAF_FOLDER_NAME/$OUTPUTNAME
	
done