# 04_PD_vs_Control_comparison.R
# Global comparison of inferred PD and Control communication networks.

source("scripts/00_setup.R")

cellchat_pd <- readRDS(file.path(OUTPUT_DIR, "cellchat_PD_full_network.rds"))
cellchat_ctrl <- readRDS(file.path(OUTPUT_DIR, "cellchat_Control_full_network.rds"))

object.list <- list(PD = cellchat_pd, Control = cellchat_ctrl)

cellchat_merged <- mergeCellChat(
  object.list,
  add.names = names(object.list)
)

saveRDS(
  cellchat_merged,
  file.path(OUTPUT_DIR, "cellchat_PD_vs_Control_merged.rds")
)

# Cell-type-pair network matrices.
pd_net <- cellchat_pd@net$weight
ctrl_net <- cellchat_ctrl@net$weight

diff_net <- pd_net - ctrl_net

saveRDS(pd_net, file.path(OUTPUT_DIR, "matrices/PD_network_weight_matrix.rds"))
saveRDS(ctrl_net, file.path(OUTPUT_DIR, "matrices/Control_network_weight_matrix.rds"))
saveRDS(diff_net, file.path(OUTPUT_DIR, "matrices/PD_minus_Control_network_matrix.rds"))

# Number of inferred interactions per cell-type pair.
pd_count <- cellchat_pd@net$count
ctrl_count <- cellchat_ctrl@net$count
count_diff <- pd_count - ctrl_count

saveRDS(pd_count, file.path(OUTPUT_DIR, "matrices/PD_network_count_matrix.rds"))
saveRDS(ctrl_count, file.path(OUTPUT_DIR, "matrices/Control_network_count_matrix.rds"))
saveRDS(count_diff, file.path(OUTPUT_DIR, "matrices/PD_minus_Control_count_matrix.rds"))

# Global network summaries used for descriptive comparison.
global_summary <- data.frame(
  condition = c("PD", "Control"),
  interaction_pairs = c(
    sum(pd_count > 0),
    sum(ctrl_count > 0)
  ),
  total_interaction_count = c(
    sum(pd_count),
    sum(ctrl_count)
  ),
  total_interaction_weight = c(
    sum(pd_net),
    sum(ctrl_net)
  )
)

write.csv(
  global_summary,
  file.path(OUTPUT_DIR, "tables/global_network_summary.csv"),
  row.names = FALSE
)

# Pathway-level comparison.
pathway_pd <- cellchat_pd@netP$pathways
pathway_ctrl <- cellchat_ctrl@netP$pathways

write.csv(
  data.frame(pathway = pathway_pd),
  file.path(OUTPUT_DIR, "tables/PD_detected_pathways.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(pathway = pathway_ctrl),
  file.path(OUTPUT_DIR, "tables/Control_detected_pathways.csv"),
  row.names = FALSE
)

# CellChat comparison plots can be generated interactively:
# compareInteractions(cellchat_merged, group = c(1, 2))
# rankNet(cellchat_merged, mode = "comparison")
