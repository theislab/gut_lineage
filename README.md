# Cell fate priming via non-canonical Wnt/PCP signalling in the small intestinal epithelium - lineage analysis in single-cell RNAsequencing data

This repository contains all scripts to reproduce the results of the single-cell data from:
A. Böttcher, M. Büttner, S. Tritschler, M. Sterr, A. Aliluev, ..., F.J.Theis, H. Lickert, "Non-canonical Wnt/PCP signalling regulates intestinal stem cell lineage priming towards enteroendocrine and Paneth cell fates", Nature Cell Biology, January 2021 - see manuscript [here](https://www.nature.com/articles/s41556-020-00617-2) and [author correction here](https://www.nature.com/articles/s41556-021-00667-0).

The notebooks contain code for the following analyses:

* QC, preprocessing -> gut_AB_AL_cell_filtering.ipynb
* Normalisation and preparation of batch effect correction -> gut_AB_AL_preBatch_cor.ipynb 
* Batch correction using adjusted ComBat -> gut_AB_AL_batch_cor.ipynb
* Cell type annotation (control samples)
  * ISC: gut_AB_AL_control_ISC_annotation.ipynb
  * Goblet, Paneth and Tuft cells: gut_AB_AL_GPT_annotation.ipynb
  * EEC: gut_AB_AL_EEC_annotation.ipynb
* Analysis of control groups and pseudotime analysis -> gut_AB_AL_cell_identity-lineage_inference.ipynb
* Annotation and analysis of mutant groups (Celsr1 crsh hemizygous mouse line) -> gut_AB_AL_mutant_analysis.ipynb

Differential expression tests were carried out using limma in R.

In addition, we added several Python functions, which we used to carry out the analysis:
* cal_density.py: Computes a density plot on an embedding, now part of `scanpy` in `tl.embedding`
* combat.py: Adjusted ComBat function that accounts for changes in cell type composition for the batch correction
* bar_frequency.py: Creates a barplot of the composition for a certain covariate across conditions
* genes_to_xls.py: Enables saving of the `tl.rank_genes_groups` results as Excel table 

The data has been deposited in GEO under accession number [GSE152325](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSE152325). The preprocessed, filtered and annotated count matrices are provided as supplementary file as a Anndata object (h5ad-file).

For further exploration load the `adata.h5ad` into a cellxgene browser for visualization or into a python-session for additional analyses using `scanpy`.

Note that the analysis was done with `scanpy v. 1.3.1`. Some functions have changed in newer versions of scanpy. For other package versions please consult the notebook or the methods in the supplementary information of the manuscript. Numeric results can vary depending on package versions and e.g. affect clustering.

If the materials in this repository are of use to you, please cite publication linked above.
