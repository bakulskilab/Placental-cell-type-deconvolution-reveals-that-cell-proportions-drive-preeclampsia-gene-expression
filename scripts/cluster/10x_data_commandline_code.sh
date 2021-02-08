# Unzip and centralize filtered count matrices in one place
gunzip -c /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_40_1/outs/filtered_feature_bc_matrix/matrix.mtx.gz \
> /scratch/bakulski_root/bakulski1/kyleac/cellranger/filtered_10x_output/kc.40.1

gunzip -c /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_40_2/outs/filtered_feature_bc_matrix/matrix.mtx.gz \
> /scratch/bakulski_root/bakulski1/kyleac/cellranger/filtered_10x_output/kc.40.2

gunzip -c /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_42_1/outs/filtered_feature_bc_matrix/matrix.mtx.gz \
> /scratch/bakulski_root/bakulski1/kyleac/cellranger/filtered_10x_output/kc.42.1

gunzip -c /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_kc_42_2/outs/filtered_feature_bc_matrix/matrix.mtx.gz \
> /scratch/bakulski_root/bakulski1/kyleac/cellranger/filtered_10x_output/kc.42.2

gunzip -c /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_478/outs/filtered_feature_bc_matrix/matrix.mtx.gz \
> /scratch/bakulski_root/bakulski1/kyleac/cellranger/filtered_10x_output/pr.478

gunzip -c /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_481/outs/filtered_feature_bc_matrix/matrix.mtx.gz \
> /scratch/bakulski_root/bakulski1/kyleac/cellranger/filtered_10x_output/pr.481

gunzip -c /scratch/bakulski_root/bakulski1/kyleac/cellranger/run_cellranger_count/batch_pique_regi_tnl_484/outs/filtered_feature_bc_matrix/matrix.mtx.gz \
> /scratch/bakulski_root/bakulski1/kyleac/cellranger/filtered_10x_output/pr.484
