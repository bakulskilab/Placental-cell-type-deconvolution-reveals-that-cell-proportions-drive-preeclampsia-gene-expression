### CibersortX Docker Powershell code ###
# Estimate cell type fractions; KCPR single-cell input and microarray GSE75010 mixture
# ` allows for multi-line entry into Powershell

# Check CibersortX version and backtick multi-line
docker `
-v

# 1) Deconvolute GSE75010
docker run `
-v G:\My` Drive\Placenta_Cell_Types\RNA\placenta_cell_types_rna\data\cibersortx_local:/src/data/ `
-v G:\My` Drive\Placenta_Cell_Types\RNA\placenta_cell_types_rna\results\cibersortx_local\gse75010:/src/outdir `
cibersortx/fractions `
--username kyleac@umich.edu `
--token $CIBERSORTX_USER_TOKEN `
--verbose TRUE `
--fraction 0 `
--single_cell TRUE `
--refsample KC_Pique-Regi_sc_signature_matrix_input_2020-12-30.txt `
--mixture GSE75010_input_mixture_2020-10-13.txt `
--rmbatchSmode TRUE

docker run `
-v G:\My` Drive\Placenta_Cell_Types\RNA\placenta_cell_types_rna\data\cibersortx_local:/src/data/ `
-v G:\My` Drive\Placenta_Cell_Types\RNA\placenta_cell_types_rna\results\cibersortx_local\gse75010_gmax2000:/src/outdir `
cibersortx/fractions `
--username kyleac@umich.edu `
--token $CIBERSORTX_USER_TOKEN `
--verbose TRUE `
--G.max 2000 `
--fraction 0 `
--single_cell TRUE `
--refsample KC_Pique-Regi_sc_signature_matrix_input_2020-12-30.txt `
--mixture GSE75010_input_mixture_2020-10-13.txt `
--rmbatchSmode TRUE


# 2) Deconvolute FACS cell fractions, CibersortX recommends fraction 0 argument for droplet-based scRNA-seq
docker run `
-v G:\My` Drive\Placenta_Cell_Types\RNA\placenta_cell_types_rna\data\cibersortx_local:/src/data/ `
-v G:\My` Drive\Placenta_Cell_Types\RNA\placenta_cell_types_rna\results\cibersortx_local\facs_mixtures:/src/outdir `
cibersortx/fractions `
--username kyleac@umich.edu `
--token 027ef63505e8202ff51c424562d3c328 `
--perm 100 `
--verbose TRUE `
--fraction 0 `
--single_cell TRUE `
--refsample KC_Pique-Regi_sc_signature_matrix_input_2020-12-30.txt `
--mixture 2021-11-17_bulk_mixture_CIBERSORTx_input_fpm.txt `
--rmbatchSmode TRUE


# Deconvolute FACS cell fractions, 2000 G.max
docker run `
-v G:\My` Drive\Placenta_Cell_Types\RNA\placenta_cell_types_rna\data\cibersortx_local:/src/data/ `
-v G:\My` Drive\Placenta_Cell_Types\RNA\placenta_cell_types_rna\results\cibersortx_local\facs_mixtures_gmax2000:/src/outdir `
cibersortx/fractions `
--username kyleac@umich.edu `
--token 027ef63505e8202ff51c424562d3c328 `
--verbose TRUE `
--perm 100 `
--fraction 0 `
--single_cell TRUE `
--G.max 2000 `
--refsample KC_Pique-Regi_sc_signature_matrix_input_2020-12-30.txt `
--mixture 2021-11-17_bulk_mixture_CIBERSORTx_input_fpm.txt `
--rmbatchSmode TRUE

