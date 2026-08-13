# 05_celltype_specific_analysis.R
# Cell-type- and pathway-specific analyses used for the thesis figures.
# No Cicero or peak/promoter annotation is performed here.

source("scripts/00_setup.R")

cellchat_pd <- readRDS(file.path(OUTPUT_DIR, "cellchat_PD_full_network.rds"))
cellchat_ctrl <- readRDS(file.path(OUTPUT_DIR, "cellchat_Control_full_network.rds"))

# Helper function for extracting communication between selected populations.
get_pair <- function(object, source, target) {
  subsetCommunication(
    object,
    sources.use = source,
    targets.use = target
  )
}

# -------------------------------------------------------------------------
# Astrocytes -> OPALIN+ oligodendrocytes
# -------------------------------------------------------------------------
pd_astro_opal <- get_pair(
  cellchat_pd,
  "Astrocytes",
  "Oligodendrocytes OPALIN+"
)

ctrl_astro_opal <- get_pair(
  cellchat_ctrl,
  "Astrocytes",
  "Oligodendrocytes OPALIN+"
)

write.csv(
  pd_astro_opal,
  file.path(OUTPUT_DIR, "tables/PD_Astrocytes_to_OPALIN.csv"),
  row.names = FALSE
)
write.csv(
  ctrl_astro_opal,
  file.path(OUTPUT_DIR, "tables/Control_Astrocytes_to_OPALIN.csv"),
  row.names = FALSE
)

astro_opal_summary <- bind_rows(
  aggregate(prob ~ pathway_name, pd_astro_opal, sum) %>%
    mutate(condition = "PD"),
  aggregate(prob ~ pathway_name, ctrl_astro_opal, sum) %>%
    mutate(condition = "Control")
)

write.csv(
  astro_opal_summary,
  file.path(OUTPUT_DIR, "tables/Astrocytes_to_OPALIN_pathway_summary.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Microglia -> OPALIN+ oligodendrocytes
# -------------------------------------------------------------------------
pd_micro_opal <- get_pair(
  cellchat_pd,
  "Microglia",
  "Oligodendrocytes OPALIN+"
)

ctrl_micro_opal <- get_pair(
  cellchat_ctrl,
  "Microglia",
  "Oligodendrocytes OPALIN+"
)

write.csv(
  pd_micro_opal,
  file.path(OUTPUT_DIR, "tables/PD_Microglia_to_OPALIN.csv"),
  row.names = FALSE
)
write.csv(
  ctrl_micro_opal,
  file.path(OUTPUT_DIR, "tables/Control_Microglia_to_OPALIN.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Interneurons -> DRD1+ and DRD2+ MSNs
# -------------------------------------------------------------------------
interneurons <- c(
  "Interneurons CHRM2+ ",
  "Interneurons LHX6+ ",
  "Interneurons VWC2+"
)

pd_drd1 <- subsetCommunication(
  cellchat_pd,
  targets.use = "MSNs DRD1+"
)
ctrl_drd1 <- subsetCommunication(
  cellchat_ctrl,
  targets.use = "MSNs DRD1+"
)

pd_drd2 <- subsetCommunication(
  cellchat_pd,
  targets.use = "MSNs DRD2+"
)
ctrl_drd2 <- subsetCommunication(
  cellchat_ctrl,
  targets.use = "MSNs DRD2+"
)

pd_drd1_inter <- subset(pd_drd1, source %in% interneurons)
ctrl_drd1_inter <- subset(ctrl_drd1, source %in% interneurons)

pd_drd2_inter <- subset(pd_drd2, source %in% interneurons)
ctrl_drd2_inter <- subset(ctrl_drd2, source %in% interneurons)

write.csv(
  pd_drd1_inter,
  file.path(OUTPUT_DIR, "tables/PD_Interneurons_to_DRD1.csv"),
  row.names = FALSE
)
write.csv(
  ctrl_drd1_inter,
  file.path(OUTPUT_DIR, "tables/Control_Interneurons_to_DRD1.csv"),
  row.names = FALSE
)
write.csv(
  pd_drd2_inter,
  file.path(OUTPUT_DIR, "tables/PD_Interneurons_to_DRD2.csv"),
  row.names = FALSE
)
write.csv(
  ctrl_drd2_inter,
  file.path(OUTPUT_DIR, "tables/Control_Interneurons_to_DRD2.csv"),
  row.names = FALSE
)

# Summarise pathway-level communication probability.
summarise_pathways <- function(pd_df, ctrl_df) {
  bind_rows(
    aggregate(prob ~ pathway_name, pd_df, sum) %>%
      mutate(condition = "PD"),
    aggregate(prob ~ pathway_name, ctrl_df, sum) %>%
      mutate(condition = "Control")
  ) %>%
    arrange(pathway_name, condition)
}

write.csv(
  summarise_pathways(pd_drd1_inter, ctrl_drd1_inter),
  file.path(OUTPUT_DIR, "tables/Interneurons_to_DRD1_pathway_summary.csv"),
  row.names = FALSE
)

write.csv(
  summarise_pathways(pd_drd2_inter, ctrl_drd2_inter),
  file.path(OUTPUT_DIR, "tables/Interneurons_to_DRD2_pathway_summary.csv"),
  row.names = FALSE
)
