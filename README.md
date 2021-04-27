# Gut lineage analysis

This repository contains all scripts to reproduce the results of the single-cell data from:
A. Böttcher, M. Büttner, S. Tritschler, M. Sterr, A. Aliluev, ..., F.J.Theis, H. Lickert, "Non-canonical Wnt/PCP signalling regulates intestinal stem cell lineage priming towards enteroendocrine and Paneth cell fates", Nature Cell Biology, January 2021 - see manuscript [here](https://www.nature.com/articles/s41556-020-00617-2) and [author correction here](https://www.nature.com/articles/s41556-021-00667-0).

The notebooks contain code for the following analyses:

* QC, preprocessing, normalisation 
* Batch correction using adjusted ComBat
* Clustering and annotation steps (input data are raw count matrices)
* Cell type annotation
* Analysis of control groups 
* Analysis of mutant groups (Celsr1 crsh hemizygous mouse line)
* Pseudotime analysis

Differential expression tests were carried out using limma in R.

The data has been deposited in GEO under accession number [GSE152325](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSE152325). The preprocessed, filtered and annotated count matrices are provided as supplementary file as a Anndata object (h5ad-file).

For further exploration load the `adata.h5ad` into a cellxgene browser for visualization or into a python-session for additional analyses using `scanpy`.

Note that the analysis was done with `scanpy v. 1.3.1`. Some functions have changed in newer versions of scanpy. For other package versions please consult the notebook or the methods in the supplementary information of the manuscript. Numeric results can vary depending on package versions and e.g. affect clustering.

If the materials in this repo are of use to you, please consider citing the above publication.
