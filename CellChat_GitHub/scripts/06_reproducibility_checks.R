# 06_reproducibility_checks.R
# Compact checks corresponding to the CellChat objects used in the thesis.

source("scripts/00_setup.R")

cellchat_pd <- readRDS(file.path(OUTPUT_DIR, "cellchat_PD_full_network.rds"))
cellchat_ctrl <- readRDS(file.path(OUTPUT_DIR, "cellchat_Control_full_network.rds"))

cat("CellChat version:", as.character(packageVersion("CellChat")), "\n")
cat("PD cells:", ncol(cellchat_pd@data.signaling), "\n")
cat("Control cells:", ncol(cellchat_ctrl@data.signaling), "\n")

cat("PD cell types:\n")
print(levels(cellchat_pd@idents))

cat("Control cell types:\n")
print(levels(cellchat_ctrl@idents))

cat("PD non-zero cell-type interactions:",
    sum(cellchat_pd@net$count > 0), "\n")
cat("Control non-zero cell-type interactions:",
    sum(cellchat_ctrl@net$count > 0), "\n")

cat("PD total interaction count:",
    sum(cellchat_pd@net$count), "\n")
cat("Control total interaction count:",
    sum(cellchat_ctrl@net$count), "\n")

cat("PD total interaction weight:",
    sum(cellchat_pd@net$weight), "\n")
cat("Control total interaction weight:",
    sum(cellchat_ctrl@net$weight), "\n")

cat("PD detected pathways:",
    length(cellchat_pd@netP$pathways), "\n")
cat("Control detected pathways:",
    length(cellchat_ctrl@netP$pathways), "\n")

cat("PD pathways:\n")
print(cellchat_pd@netP$pathways)

cat("Control pathways:\n")
print(cellchat_ctrl@netP$pathways)
