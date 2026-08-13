# CellChat analysis of PD vs Control

This repository contains the CellChat analysis used to investigate inferred cell-cell communication changes between Parkinson's disease (PD) and control conditions from single-cell H3K27ac CUT&Tag-derived gene activity profiles.

## Scope

This repository contains **CellChat analysis only**.

It does not contain:
- Cicero analysis
- peak-to-gene or peak/promoter annotation
- raw sequencing data
- the complete Seurat object
- HPC/R package libraries

## Analysis overview

1. Load the preprocessed Seurat object containing the `gene_activity` assay.
2. Split cells into PD and control conditions.
3. Extract the `gene_activity` assay (`layer = "data"`).
4. Construct separate CellChat objects for PD and control.
5. Use the human CellChatDB.
6. Identify overexpressed signaling genes and ligand-receptor interactions.
7. Infer communication probabilities.
8. Filter inferred communication.
9. Aggregate interactions into signaling pathways and cell-cell networks.
10. Compute pathway-level centrality.
11. Compare PD and control communication networks.
12. Perform selected cell-type-specific analyses.

## Reproducibility

The verified CellChat environment used for the thesis included:

- R 4.4.3
- CellChat 1.6.1
- Seurat 5.4.0
- Signac 1.16.0
- Matrix 1.7.5
- ggplot2 4.0.2
- dplyr 1.2.0

CellChat inference parameters verified from the saved objects:

- `type = "triMean"`
- `trim = 0.1`
- `raw.use = TRUE`
- `population.size = FALSE`
- `nboot = 100`
- `seed.use = 1`
- `Kh = 0.5`
- `n = 1`

The CellChat objects used in the thesis contained 14 annotated cell populations.

## Important note about input data

The raw/preprocessed Seurat object is not included in this repository. Before running the scripts, edit `INPUT_RDS` in `scripts/00_setup.R` to point to the local copy of the appropriate Seurat object.

The input object must contain:
- a `gene_activity` assay
- `pd_status` metadata
- `celltype` metadata

## Running the analysis

Run the scripts in order:

```text
scripts/00_setup.R
scripts/01_prepare_input.R
scripts/02_create_cellchat.R
scripts/03_infer_communication.R
scripts/04_PD_vs_Control_comparison.R
scripts/05_celltype_specific_analysis.R
scripts/06_reproducibility_checks.R
```

## Main cell-type-specific analyses

The repository includes code for:

- Astrocytes -> Oligodendrocytes OPALIN+
- Microglia -> Oligodendrocytes OPALIN+
- Interneurons -> MSNs DRD1+
- Interneurons -> MSNs DRD2+

## Interpretation

CellChat results represent **inferred** ligand-receptor-mediated communication based on the input molecular profiles and the CellChatDB prior-knowledge database. They should not be interpreted as direct experimental measurements of physical cell-cell communication.

## License

Add the licence required by the project/supervisors before publication.
