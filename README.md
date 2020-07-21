# gut_lineage

This repository contains all scripts to reproduce the results of the single-cell data from:
A. Böttcher et al., "Wnt/PCP-primed intestinal stem cells differentiate into enteroendocrine or Paneth cells", 2020

The notebooks contain code for the following analyses:

* QC, preprocessing, clustering and annotation steps (input data are raw count matrices)
* Data annotation
* Analysis of control groups
* Pseudotime analysis

Differential expression tests were carried out using limma in R.

The data has been deposited in GEO under accession number GSE152325. The preprocessed, filtered and annotated count matrices are provided as supplementary file as a Anndata object (h5ad-file).

For further exploration load the adata.h5ad into a cellxgene browser for visualization or into a python-session for additional analyses using scanpy.

Note that the analysis was done with scanpy v. 1.3.1. Some functions have changed in newer versions of scanpy. For other package versions please consult the notebook or the methods in the supplementary information of the manuscript. Numeric results can vary depending on package versions and e.g. affect clustering.

If the materials in this repo are of use to you, please consider citing the above publication.
