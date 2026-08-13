# 01_prepare_input.R
# Prepare gene-activity matrices and metadata for CellChat.

source("scripts/00_setup.R")

new_obj <- readRDS(INPUT_RDS)

if (!"gene_activity" %in% names(new_obj@assays)) {
  stop("The Seurat object does not contain a 'gene_activity' assay.")
}
if (!all(c("pd_status", "celltype") %in% colnames(new_obj[[]]))) {
  stop("Metadata must contain 'pd_status' and 'celltype'.")
}

DefaultAssay(new_obj) <- "gene_activity"

obj_pd <- subset(new_obj, subset = pd_status == "pd")
obj_ctrl <- subset(new_obj, subset = pd_status == "control")

ga_pd <- GetAssayData(obj_pd, assay = "gene_activity", layer = "data")
ga_ctrl <- GetAssayData(obj_ctrl, assay = "gene_activity", layer = "data")

meta_pd <- obj_pd[[]]
meta_ctrl <- obj_ctrl[[]]

# Keep only cells with the expected cell-type labels.
keep_pd <- meta_pd$celltype %in% CELLTYPES
keep_ctrl <- meta_ctrl$celltype %in% CELLTYPES

ga_pd <- ga_pd[, keep_pd, drop = FALSE]
ga_ctrl <- ga_ctrl[, keep_ctrl, drop = FALSE]
meta_pd <- meta_pd[keep_pd, , drop = FALSE]
meta_ctrl <- meta_ctrl[keep_ctrl, , drop = FALSE]

meta_pd <- meta_pd[colnames(ga_pd), , drop = FALSE]
meta_ctrl <- meta_ctrl[colnames(ga_ctrl), , drop = FALSE]

stopifnot(identical(rownames(meta_pd), colnames(ga_pd)))
stopifnot(identical(rownames(meta_ctrl), colnames(ga_ctrl)))

saveRDS(ga_pd, file.path(OUTPUT_DIR, "ga_pd_gene_activity.rds"))
saveRDS(ga_ctrl, file.path(OUTPUT_DIR, "ga_ctrl_gene_activity.rds"))
saveRDS(meta_pd, file.path(OUTPUT_DIR, "meta_pd.rds"))
saveRDS(meta_ctrl, file.path(OUTPUT_DIR, "meta_ctrl.rds"))

message("PD: ", ncol(ga_pd), " cells; Control: ", ncol(ga_ctrl), " cells")
