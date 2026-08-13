# CellChat analysis of PD vs Control
# 00_setup.R
# Purpose: load packages, define reproducibility settings and output paths.

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
})

set.seed(1)

# Edit this path for your local/HPC environment.
INPUT_RDS <- "/scratch/prj/bcn_marzi_lab/analysis_cutandtag_pd_sc/student_data_package/data_in/singlecell_data/merged_clustered_sc_data_clean_YIFEI.rds "

OUTPUT_DIR <- "results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "matrices"), showWarnings = FALSE, recursive = TRUE)

# Cell type labels used in the thesis.
CELLTYPES <- c(
  "Astrocytes",
  "Interneurons CHRM2+ ",
  "Interneurons LHX6+ ",
  "Interneurons VWC2+",
  "Microglia",
  "MSNs DRD1+",
  "MSNs DRD2+",
  "MSNs GRIK3+",
  "MSNs NPNT+ ",
  "MSNs RXFP1+",
  "Oligodendrocyte precursor cells (OPCs)",
  "Oligodendrocytes OPALIN+",
  "Oligodendrocytes PLECHG1+",
  "Perivascular fibroblast-like"
)

stopifnot(packageVersion("CellChat") == "1.6.1")
