# 03_infer_communication.R
# CellChat ligand-receptor and pathway inference.
#
# Parameters reproduce the verified CellChat objects used in the thesis:
# type.mean = triMean; trim = 0.1; raw.use = TRUE;
# population.size = FALSE; nboot = 100; seed.use = 1;
# Kh = 0.5; n = 1.

source("scripts/00_setup.R")

cellchat_pd <- readRDS(file.path(OUTPUT_DIR, "cellchat_PD_initial.rds"))
cellchat_ctrl <- readRDS(file.path(OUTPUT_DIR, "cellchat_Control_initial.rds"))

run_cellchat <- function(x) {
  x <- subsetData(x)

  x <- identifyOverExpressedGenes(x)
  x <- identifyOverExpressedInteractions(x)

  x <- computeCommunProb(
    x,
    type = "triMean",
    trim = 0.1,
    raw.use = TRUE,
    population.size = FALSE,
    nboot = 100,
    seed.use = 1,
    Kh = 0.5,
    n = 1
  )

  # Retain communication supported by at least 10 cells.
  x <- filterCommunication(x, min.cells = 10)

  x <- computeCommunProbPathway(x)

  x <- aggregateNet(x)

  x <- netAnalysis_computeCentrality(x, slot.name = "netP")

  x
}

cellchat_pd <- run_cellchat(cellchat_pd)
cellchat_ctrl <- run_cellchat(cellchat_ctrl)

saveRDS(
  cellchat_pd,
  file.path(OUTPUT_DIR, "cellchat_PD_full_network.rds")
)
saveRDS(
  cellchat_ctrl,
  file.path(OUTPUT_DIR, "cellchat_Control_full_network.rds")
)
