# 02_create_cellchat.R
# Create PD and Control CellChat objects and assign the human CellChatDB.

source("scripts/00_setup.R")

ga_pd <- readRDS(file.path(OUTPUT_DIR, "ga_pd_gene_activity.rds"))
ga_ctrl <- readRDS(file.path(OUTPUT_DIR, "ga_ctrl_gene_activity.rds"))
meta_pd <- readRDS(file.path(OUTPUT_DIR, "meta_pd.rds"))
meta_ctrl <- readRDS(file.path(OUTPUT_DIR, "meta_ctrl.rds"))

CellChatDB <- CellChatDB.human

cellchat_pd <- createCellChat(
  object = ga_pd,
  meta = meta_pd,
  group.by = "celltype"
)
cellchat_ctrl <- createCellChat(
  object = ga_ctrl,
  meta = meta_ctrl,
  group.by = "celltype"
)

cellchat_pd@DB <- CellChatDB
cellchat_ctrl@DB <- CellChatDB

saveRDS(cellchat_pd, file.path(OUTPUT_DIR, "cellchat_PD_initial.rds"))
saveRDS(cellchat_ctrl, file.path(OUTPUT_DIR, "cellchat_Control_initial.rds"))
