library(Seurat)
library(SeuratDisk)
library(CellChat)
axolotl <- readRDS("G:/DPayz_sysAct_seurat_20230327.rds")
axolottest <- axolotl
axolotl_contra <- subset(axolottest, subset = sample == c('contra'))
axolotl_intact <- subset(axolottest, subset = sample == c('intact'))
head(axolotl_intact)
table(axolotl_intact@active.ident)

interaction_input <- read.csv(file = 'C://Users/Emil/10X/cellchat/interaction_input_CellChatDB.csv', row.names = 1)
complex_input <- read.csv(file = 'C://Users/Emil/10X/cellchat/complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = 'C://Users/Emil/10X/cellchat/cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = 'C://Users/Emil/10X/cellchat/geneInfo_input_CellChatDB.csv', row.names = 1)

CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo
CellChatDB.mouse <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <-IdentifyOverExpressedGenes(cellchat)
cellchat@var.features
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
