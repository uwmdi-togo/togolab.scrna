# togolab.scrna

Single-cell RNA-seq and downstream -omics analysis helpers for the
Bjornstad-Pyle-Tommerdahl (Togo) lab. Companion to the lightweight
[`togolab`](https://github.com/uwmdi-togo/togolab) core package — the heavy
analysis dependencies (Seurat, nebula, fgsea, slingshot, …) live here so they
don't burden the core.

## Install

```r
# install.packages("remotes")
remotes::install_github("uwmdi-togo/togolab.scrna")
library(togolab.scrna)
```

This pulls in `togolab` automatically (declared in `Remotes:`). Each function
only needs its specific heavy package at call time, and errors with a clear
install message if it's missing.

## Functions

These work for both the ATTEMPT (treatment x visit) and PB90 (disease-group)
designs — analysis is formula-/column-driven, with no design hardcoded.

NEBULA differential expression: `togo_run_nebula()` (single formula-driven fit),
`togo_run_nebula_parallel()` (per-gene in parallel), `togo_process_nebula_results()`.

scRNA wrangling/QC: `togo_make_subsets()`, `togo_run_doubletfinder()`.

Feature UMAPs & composition: `togo_prepare_umap_metadata()`,
`togo_plot_feature_umap()`, `togo_celltype_proportions()`, `togo_celltype_pie()`.

Volcano & pathways: `togo_plot_volcano()`, `togo_prepare_gmt()`,
`togo_matrix_to_list()`, `togo_clean_pathway_names()`,
`togo_filter_redundant_pathways()`, `togo_plot_fgsea()`.

Trajectory (Slingshot): `togo_slingshot_setup()`, `togo_run_slingshot()`.

## Example

```r
library(togolab.scrna)

# Differential expression for a disease-group contrast (PB90):
fit <- togo_run_nebula(~ group, object = so, id_col = "record_id")

# ...or a treatment x visit interaction (ATTEMPT):
res <- togo_run_nebula_parallel(so, formula = ~ visit * treatment,
                                id_col = "subject_id")
summ <- togo_process_nebula_results(res$results)
```

## Notes

Not included (too ATTEMPT-specific to generalize): the SomaScan/limma proteomics
functions and the PRE/POST-by-treatment pseudotime *plots*. Ask if you want
generalized versions.

After adding/editing functions, run `devtools::document()` to regenerate the
`man/` help pages and `NAMESPACE` from the roxygen comments.
