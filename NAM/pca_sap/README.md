This data was produced by running the following script, a wrapper for pca.R:
/pl/active/Morris_CSU/Ted_Monyak/NAM/run_pca.sh

Genotype data came from:
/pl/active/Morris_CSU/GBS_data/NAM_Hu_Faye_Lasky_GBS_V5_updated/NAM/WGS/Ted_NAM_SAP/NAM_SAP_GBS_shared_chrid_refnorm.vcf.gz

PCA was run on the SAP data only
The NAM, plus the founder lines, were projected onto the SAP

Files:
full_pca_projection.rds: All of the PCs of the NAM projections
pca_df.rds: A dataframe containing the NAM projections onto the first 9 PCs
vaf.rds: The variance explained by each PC in the SAP