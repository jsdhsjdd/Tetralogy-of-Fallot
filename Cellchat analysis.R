library(patchwork)
options(stringsAsFactors = FALSE)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(devtools)
library(CellChat)
Other1 <- BGPC[,BGPC@meta.data[["celltype"]] %in% c("high_NRG1",	"GPC","OL",	"high_MKI67")]
Other2 <- ExN_data[,ExN_data@meta.data[["celltype"]] %in% c("ExN_TLE4", "ExN_CNTNAP2",	"ExN_ENC1",	"ExN_CNTN5",	"ExN_TUBB3",	"ExN_MEIS2", "ExN_POU6F2")]
Other3 <- InN_data[,InN_data@meta.data[["celltype"]] %in% c("InN_SOX4","InN_ADARB2",	"InN_RELN",	"InN_SGCZ",	"InN_MAF", "InN_CCK")]
Other4 <- GC_data[,GC_data@meta.data[["celltype"]] %in% c("Endo", "Peri",'MG')]
scRNA <- merge(Other1,c(Other2,Other3,Other4))
save(scRNA,file="scRNA.RData")
obj.list<-SplitObject(scRNA, split.by = "SampleID")
##TOF
data.input  <- obj.list[["TOF"]]@assays[["RNA"]]@data
celltype  <- obj.list[["TOF"]]@meta.data[["celltype"]]
identity = data.frame(group = obj.list[["TOF"]]$celltype, row.names = names(obj.list[["TOF"]]$celltype))
names(identity)
head(identity)
unique(identity$group)
table(identity$group)
meta <- data.frame(labels = obj.list[["TOF"]]$celltype, row.names = names(identify))
cellchat <- createCellChat(object = data.input,meta = meta, group.by = "labels")
cellchat
summary(cellchat)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
head(cellchat@meta)
cellchat <- setIdent(cellchat, ident.use = "labels") 
CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat,type = "triMean", population.size = F)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat,"SS_TOF.cellchat.rds")
df.net <- subsetCommunication(cellchat)
write.csv(df.net,file = "SS_df.net_TOF.csv",row.names = F)
groupSize <- as.numeric(table(SS_TOF.cellchat@idents))
pdf("SS_TOF_cellchat_numberandstrength.pdf", width = 10,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(SS_TOF.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(SS_TOF.cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- SS_TOF.cellchat@net$weight
pdf("SS_TOF_single_celltype.pdf", width = 20,height = 25)
par(mfrow = c(5,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

##Con.21W
data.input <- obj.list[["Con_21W"]]@assays[["RNA"]]@data
celltype  <- obj.list[["Con_21W"]]@meta.data[["celltype"]]
identity = data.frame(group = obj.list[["Con_21W"]]$celltype, row.names = names(obj.list[["Con_21W"]]$celltype))
head(identity)
unique(identity$group)
table(identity$group)
meta <- data.frame(labels = obj.list[["Con_21W"]]$celltype, row.names = names(identify))
cellchat <- createCellChat(object = data.input,meta = meta, group.by = "labels")
cellchat
summary(cellchat)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")

cellchat <- setIdent(cellchat, ident.use = "labels") 
CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat,type = "triMean", population.size = F)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat,"SS_Con_21W.cellchat.rds")
df.net <- subsetCommunication(cellchat)
write.csv(df.net,file = "SS_df.net_Con_21W.csv",row.names = F)
groupSize <- as.numeric(table(SS_Con_21W.cellchat@idents))
pdf("SS_Con.21W_numberandstrength.pdf", width = 10,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(SS_Con_21W.cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(SS_Con_21W.cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- SS_Con_21W.cellchat@net$weight
pdf("SS_Con.21W_single_celltype.pdf", width = 20,height = 25)
par(mfrow = c(5,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

#Comparative analysis
obj.list <- list(Con_21W = SS_Con_21W.cellchat,TOF =SS_TOF.cellchat)
names (obj.list)
cellchat <- mergeCellChat(obj.list,add.names = names(obj.list),cell.prefix = T)
saveRDS(cellchat, file="SS_cellchat.rds")
gg1 <- compareInteractions(cellchat, show.legend = F, color.use = c("#80B9C8","#F2A49A"),group = c(1,2),measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, color.use = c("#80B9C8","#F2A49A"),group = c(1,2),measure = "weight")
p <- gg1+gg2
ggsave("SS_overview_number_strength.pdf", p, width = 6,height = 4)
gg1 <- rankNet(cellchat, mode = "comparison", color.use = c("#80B9C8","#F2A49A"), comparison = c(1,2),stacked = T, do.stat = T)
gg2 <- rankNet(cellchat, mode = "comparison", color.use = c("#80B9C8","#F2A49A"), comparison = c(1,2),stacked = F, do.stat = T)
p <- gg1 + gg2
ggsave("SS_compare_pathway_strength.pdf", p,width = 10,height = 5.5)
pdf("SS_Compare_signal_net.pdf", width = 10,height = 10)
par(mfrow = c(2,2))
h1 <- netVisual_diffInteraction(cellchat,edge.weight.max =15,comparison = c(1, 2), weight.scale = T,color.edge = c("#F2A49A","#80B9C8"))
h2 <- netVisual_diffInteraction(cellchat, measure = "weight",comparison = c(1, 2),vertex.size.max = 15, weight.scale = T,color.edge = c("#F2A49A","#80B9C8"))
dev.off()
par(mfrow = c(1,2))
h1 <- netVisual_heatmap(cellchat ,comparison = c(1, 2))
h2 <- netVisual_heatmap(cellchat, comparison = c(1, 2), measure = "weight")
pdf("SS_Diff_number_strength_heatmap.pdf", width = 10,height = 5.5)
h1+h2
dev.off()
library(ComplexHeatmap)
pathway.union <- union(obj.list[[1]]@netP$pathways,  obj.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(obj.list[[1]], pattern = "all", signaling = pathway.union,title = names(obj.list)[1], width = 8, height = 8)
ht2 = netAnalysis_signalingRole_heatmap(obj.list[[2]], pattern = "all", signaling = pathway.union,title = names(obj.list)[2], width = 8, height = 8)
pdf("SS_compare_signal_pattern_all.pdf", width = 16,height = 8)
draw(ht1 + ht2,ht_gap = unit(0.5, "cm"))
garbage <- dev.off()
ht4 = netAnalysis_signalingRole_heatmap(obj.list[[1]], pattern = "outgoing", signaling = pathway.union,title = names(obj.list)[1], width = 8, height = 10)
ht5 = netAnalysis_signalingRole_heatmap(obj.list[[2]], pattern = "outgoing", signaling = pathway.union,title = names(obj.list)[2], width = 8, height = 10)
pdf("SS_compare_signal_pattern_outgoing.pdf", width = 16,height = 8)
draw(ht4 + ht5 ,ht_gap = unit(0.5, "cm"))
garbage <- dev.off()
ht7 = netAnalysis_signalingRole_heatmap(obj.list[[1]], pattern = "incoming", signaling = pathway.union,title = names(obj.list)[1], width = 8, height = 10)
ht8 = netAnalysis_signalingRole_heatmap(obj.list[[2]], pattern = "incoming", signaling = pathway.union,title = names(obj.list)[2], width = 8, height = 10)
pdf("SS_compare_signal_pattern_incoming.pdf", width = 16,height = 8)
draw(ht7 + ht8 ,ht_gap = unit(0.5, "cm"))
garbage <- dev.off()
num.link <- sapply(obj.list, function(x) {rowSums(x@net$count) +
    colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
pdf("bubbleplot.pdf", width = 10,height = 5)
for (i in 1:length(obj.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(obj.list[[i]],
                                               title = names(obj.list)[i], 
                                               weight.MinMax = weight.MinMax)}
patchwork::wrap_plots(plots = gg)
dev.off()

##special pathway
#NRG pathway network
pathways.show <-c ("NRG")
weight.max <- getMaxWeight(obj.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
pdf("NRG_circle.pdf", width = 10,height = 5.5)
for (i in 1:length(obj.list)) {
  netVisual_aggregate(obj.list[[i]],signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      vertex.weight = as.numeric(table(obj.list[[i]]@idents)),
                      signaling.name = paste(pathways.show, names(obj.list)[i]))  
}
dev.off()
pdf("NRG.Con.21W.NRG1cluster_chord.pdf", width = 5,height = 5)
pathways.show <-c ("NRG")
netVisual_aggregate(SS_Con_21W.cellchat, signaling = pathways.show, layout = "chord",sources.use = c(11), targets.use = c(19,17,15,16,12,13,14))
dev.off()
cellchat1 <- SS_TOF.cellchat
pdf("NRG_contribution_tof.pdf",width=5,height = 4)
netAnalysis_contribution(cellchat1,signaling = pathways.show)
dev.off()
cellchat2 <- SS_Con_21W.cellchat
pdf("NRG_contribution_Con_21W.pdf",width=5,height = 4)
netAnalysis_contribution(cellchat2,signaling = pathways.show)
dev.off()
pairLR.NRG <- extractEnrichedLR(cellchat1, signaling = pathways.show,
                                geneLR.return = FALSE)#提取显著作用的配受体对
LR.show <- pairLR.NRG[1,] 
pdf("NRG1_ERBB3_circle_CON_21W.pdf",width=10,height = 5.5)
p1<- netVisual_individual(cellchat, signaling = pathways.show, 
                          pairLR.use = LR.show, layout = "circle")
dev.off()

p1 <- plotGeneExpression(SS_TOF.cellchat,signaling='NRG',type='violin')
ggsave("NRG_geneexpression_TOF.pdf",p1,width = 10,height=5.5)
p2 <- plotGeneExpression(SS_Con_21W.cellchat,signaling='NRG',type='violin')
ggsave("NRG_geneexpression_Con21W.pdf",p2,width = 10,height=5.5)

##PTN pathway network
pathways.show <-c ("PTN")
weight.max <- getMaxWeight(obj.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
pdf("PTN_circle.pdf", width = 10,height = 5.5)
for (i in 1:length(obj.list)) {
  netVisual_aggregate(obj.list[[i]],signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      vertex.weight = as.numeric(table(obj.list[[i]]@idents)),
                      signaling.name = paste(pathways.show, names(obj.list)[i]))  
}
dev.off()



#PDGF pathway network
pathways.show <-c ("PDGF")
weight.max <- getMaxWeight(obj.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
pdf("PDGF_chord.pdf", width = 10,height = 5.5)
for (i in 1:length(obj.list)) {
  netVisual_aggregate(obj.list[[i]],signaling = pathways.show, layout = "chord",
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      vertex.weight = as.numeric(table(obj.list[[i]]@idents)),
                      signaling.name = paste(pathways.show, names(obj.list)[i]))  
}
dev.off()

#PDGF信号通路的配受体气泡图
levels(cellchat@idents$joint)
par(mfrow = c(1,2), xpd=TRUE)
pairLR.use <- extractEnrichedLR(cellchat, signaling = "PDGF")

p <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(19,20), comparison = c(1,2), pairLR.use = pairLR.use, angle.x = 45)
p
ggsave("PDGF_LR_bubble.pdf", p, width = 4, height = 2)


p1 <- plotGeneExpression(SS_TOF.cellchat,signaling='PDGF',type='violin')
ggsave("PDGF_geneexpression_TOF.pdf",p1,width = 10,height=5.5)
p2 <- plotGeneExpression(SS_Con_21W.cellchat,signaling='PDGF',type='violin')
ggsave("PDGF_geneexpression_Con21W.pdf",p2,width = 10,height=5.5)


#VEGF pathway network
pathways.show <-c ("VEGF")
weight.max <- getMaxWeight(obj.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
pdf("VEGF_chord.pdf", width = 10,height = 5.5)
for (i in 1:length(obj.list)) {
  netVisual_aggregate(obj.list[[i]],signaling = pathways.show, layout = "chord",
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      vertex.weight = as.numeric(table(obj.list[[i]]@idents)),
                      signaling.name = paste(pathways.show, names(obj.list)[i]))  
}
dev.off()

#GAS pathway network
pathways.show <-c ("GAS")
weight.max <- getMaxWeight(obj.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
pdf("GAS_chord.pdf", width = 10,height = 5.5)
for (i in 1:length(obj.list)) {
  netVisual_aggregate(obj.list[[i]],signaling = pathways.show, layout = "chord",
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      vertex.weight = as.numeric(table(obj.list[[i]]@idents)),
                      signaling.name = paste(pathways.show, names(obj.list)[i]))  
}
dev.off()


#TGFb pathway network
pathways.show <- c("TGFb") 
pdf("TGFb.Con.21W.chord.pdf", width = 5,height = 5.5)
netVisual_aggregate(SS_Con_21W.cellchat, signaling = pathways.show, layout = "chord")
dev.off()


#COMPLEMENT pathway network
pathways.show <- c("COMPLEMENT") 
pdf("COMPLEMENT.Con.21W.chord.pdf", width = 5,height = 5.5)
netVisual_aggregate(SS_Con_21W.cellchat, signaling = pathways.show, layout = "chord")
dev.off()

