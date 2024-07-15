##subtype identification
load("obj.combined2.RData")
InN_data <- obj.combined2[,obj.combined2@meta.data[["celltype"]] %in% c("InN")]
DefaultAssay(InN_data) <- "RNA" 
colnames(InN_data@meta.data)
InN_data[['seurat_clusters']] <- NULL
InN_data[['integrated_snn_res.0.6']] <- NULL
InN_data[['celltype']] <- NULL

InN_data <-  NormalizeData(InN_data,normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 2000)
InN_data <- ScaleData(InN_data, features = rownames(InN_data))
InN_data <- RunPCA(InN_data, verbose = FALSE)
plot2 <- ElbowPlot(InN_data,ndims = 50)
dimNums = 8
InN_data <- FindNeighbors(InN_data, reduction = "pca", dims = 1:dimNums)
InN_data <- FindClusters(InN_data,resolution = 0.5)
InN_data <- RunUMAP(InN_data, dims = 1:dimNums)
InN_data <- InN_data[,InN_data@meta.data[["seurat_clusters"]] %in% c('0',"1","3", '4',"5",'6','7')]
DefaultAssay(InN_data) <- "RNA" 
colnames(InN_data@meta.data)
InN_data[['seurat_clusters']] <- NULL
InN_data[['RNA_snn_res.0.5']] <- NULL

InN_data <-  NormalizeData(InN_data,normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 2000)
InN_data <- ScaleData(InN_data, features = rownames(InN_data))
InN_data <- RunPCA(InN_data, verbose = FALSE)
plot2 <- ElbowPlot(InN_data,ndims = 50)
dimNums = 5 
InN_data <- FindNeighbors(InN_data, reduction = "pca", dims = 1:dimNums)
InN_data <- FindClusters(InN_data,resolution = 0.5)
InN_data <- RunUMAP(InN_data, dims = 1:dimNums)
markers <- FindAllMarkers(InN_data, logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
current.cluster.ids <- c("0",	"1",	"2", "3",	"4",	"5",	"6")
new.cluster.ids <- c("InN_SOX4","InN_RUNX1",	"InN_ADARB2",	"InN_RELN",	"InN_SGCZ",	"InN_MAF", "InN_CCK") 
InN_data@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(InN_data@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(InN_data@meta.data$celltype)
current.cluster.ids <- c("0",	"1",	"2", "3",	"4",	"5",	"6")
names(new.cluster.ids) <- levels(InN_data)
InN_data<- RenameIdents(InN_data, new.cluster.ids)
DimPlot(InN_data, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
color2=c(InN_SOX4= "#E5A84B", InN_RUNX1 = "#509296",InN_ADARB2 = "#EE7959",InN_RELN = "#9BB496",
         InN_MAF = "#F08080",InN_CCK = "#A64036")
p9 <- DimPlot(InN_data, cols = color2, reduction = "umap", group.by = "ident", pt.size=0.1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=1),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "UMAP.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')

library(pheatmap)
DefaultAssay(InN_data) <- "RNA"
markers <- FindAllMarkers(InN_data, logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
sig_markers <- markers %>% group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
head(sig_markers)
genes<- unique(sig_markers$gene)
aver_dt <- AverageExpression(InN_data,features = genes,group.by ='celltype',slot="data")
aver_dt<-as.data.frame(aver_dt$RNA)
aver_dt[1:6,1:6]
cell_anno<-data.frame(cell_anno=colnames(aver_dt),row.names =colnames(aver_dt))
color2=c(InN_SOX4= "#E5A84B", InN_SGCZ = "#509296",InN_ADARB2 = "#EE7959",InN_RELN = "#9BB496",
         InN_MAF = "#F08080",InN_CCK = "#A64036")
celltype_col<-color2
names(celltype_col)<-cell_anno$cell_anno
anno_col<-list(cell_anno=celltype_col)
anno_col

p1=pheatmap(as.matrix(aver_dt),scale='row',cellwidth = 10, cellheight = 5,cluster_rows = FALSE,cluster_cols =FALSE,annotation_col = cell_anno,annotation_colors = anno_col,border_color = 'white')
ggsave(filename = "heatmap+top10.pdf", plot = p1, device = 'pdf', width = 10, height = 18, units = 'cm')

InN_data$SampleID <- "SampleID"
InN_data$SampleID[ as.character(InN_data$orig.ident) %in% c("TOF") ] <- "TOF"
InN_data$SampleID[ as.character(InN_data$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
InN_data$SampleID[ as.character(InN_data$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
InN_data$SampleID[ as.character(InN_data$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
InN_data$SampleID[ as.character(InN_data$orig.ident) %in% c("Con_4") ] <- "Con.20W"
InN_data$SampleID[ as.character(InN_data$orig.ident) %in% c("Con_8")] <- "Con.23W"
InN_data$SampleID[ as.character(InN_data$orig.ident) %in% c("Con_9") ] <- "Con.25W"

p9 <- DimPlot(InN_data, cols = c("#8DD3C7", "#FDB462", "#BEBADA", "#FB8072", "#80B1D3","#FFFFB3", "#B3DE69"), reduction = "umap", group.by = "SampleID",pt.size=0.1, label = F,label.size=5, repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "InN_data+umap+SampleID.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')

AB <- InN_data[,InN_data@meta.data[["SampleID"]] %in% c("TOF")]
p9 <- DimPlot(AB, cols = color2, reduction = "umap", group.by = "celltype",pt.size=0.1, label = F,label.size=5, repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "InN_data+umap+TOF.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')

AB <- InN_data[,InN_data@meta.data[["SampleID"]] %in% c("Con.21W")]
p9 <- DimPlot(AB, cols = color2, reduction = "umap", group.by = "celltype",pt.size=0.1, label = F,label.size=5, repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "InN_data+umap+Con.21W.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')


AB <- InN_data[,InN_data@meta.data[["SampleID"]] %in% c("TOF","Con.21W")]
p9 <- DimPlot(AB, cols = c('#EAA096','#7CB2C0'), reduction = "umap", group.by = "SampleID", pt.size=0.5, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "InN_data+UMAP+TOF_Con21.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')

##Mfuzz analysis
library(Mfuzz)
library(Seurat)
library(progeny)
library(tidyr)
library(tibble)
library(dplyr)
library(viridisLite)
library(ggplot2)
library("BiocParallel")
register(SnowParam(6))
Sys.setenv("VROOM_CONNECTION_SIZE"=88888888)
scRNA<- InN_data
rm(InN_data)
scRNA$SampleID <- "SampleID"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("TOF") ] <- "TOF"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
scRNA$SampleID[ as.character(scRNA$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_4") ] <- "Con.20W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_8")] <- "Con.23W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_9") ] <- "Con.25W"
scRNA@meta.data[["SampleID"]]<-factor(scRNA@meta.data[["SampleID"]], levels=c("TOF","Con.17W","Con.20W", "Con.21W","Con.22W","Con.23W","Con.25W"))
scRNA=ScaleData(scRNA)
age.averages <- AverageExpression(scRNA,group.by = "SampleID")
df1 <- age.averages[["RNA"]]
mat <- as.matrix(df1)
dt <- new("ExpressionSet",exprs = mat)
dt.r <- filter.NA(dt, thres=0.25)
dt.f <- fill.NA(dt.r,mode="mean")
tmp <- filter.std(dt.f,min.std=0.1)
dt.s <- standardise(tmp)
df.s <- dt.s@assayData$exprs
m1 <- mestimate(dt.s)
set.seed(007)
cl <- mfuzz(dt.s,c=7,m=m1)
mfuzz.plot(dt.s,cl,mfrow=c(3,3), new.window= FALSE, time.labels=colnames(dt.s))
library(RColorBrewer)
mycol <- c("#E5A84B","#EE7959","#8EC9CF")
mycolor <- colorRampPalette(mycol)(20)
pdf("RNA_TOF_InN expression changes trends_0.1.pdf", width = 14,height = 7)
mfuzz.plot(dt.s,cl,mfrow=c(2,3), new.window= FALSE, time.labels=colnames(dt.s), colo = mycolor)
dev.off()
gene_cluster <- data.frame(df.s,cluster=cl$cluster)
write.csv(gene_cluster, 'RNA_TOF_InN expression changes trends_0.1.csv')
#Mfuzz:cluster4
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix)
library(readxl)
library(ggpubr)
BiocManager::install("ggpubr")
retina_gene<-read_xlsx("cluster4.xlsx")
gene<-as.list(retina_gene)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "retina")
colnames(AB@meta.data)[8]<-"retina_Score"
AB@meta.data[["celltype"]]<-factor(AB@meta.data[["celltype"]], levels=c("InN_SOX4","InN_ADARB2",	"InN_RELN",	"InN_SGCZ",	"InN_MAF", "InN_CCK"))

color2=c(InN_SOX4= "#E5A84B", InN_SGCZ = "#509296",InN_ADARB2 = "#EE7959",InN_RELN = "#9BB496",
         InN_MAF = "#F08080",InN_CCK = "#A64036")
library(ggpubr)
p.AddModuleScore <- ggviolin(AB@meta.data, x = "celltype", y = "retina_Score",
                             color = "celltype",add = 'mean_sd',fill = 'celltype',
                             add.params = list(color = "black")) + 
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  scale_color_manual(values = color2) + 
  scale_fill_manual(values = color2) +
  #theme(axis.text.x.bottom = element_text(angle = 0,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')
p.AddModuleScore
ggsave('RNA_TOF_InN expression changes trends_0.1_cluster4_vlo.pdf',width = 8,height = 4)


##cytotrace
devtools::install_local("D:/科研/2023年/神经管畸形/NTD/参考/CytoTRACE_0.3.3.tar.gz")
conda_create("cytoTRACE",python_version = '3.12.2')
BiocManager::install("sva")
BiocManager::install("genefilter",force = TRUE)
library(reticulate)
library(CytoTRACE)
reticulate::py_install(packages = c (  "scanoramaCT" )) 

scRNA <- InN_data
scRNA$SampleID <- "SampleID"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("TOF") ] <- "TOF"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
scRNA$SampleID[ as.character(scRNA$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_4") ] <- "Con.20W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_8")] <- "Con.23W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_9") ] <- "Con.25W"

InN.21W <- scRNA[,scRNA@meta.data[["SampleID"]] %in% c("TOF","Con.21W")]
exp1 <- as.matrix(InN.21W@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
results <- CytoTRACE(exp1,ncores = 1)
phenot <- InN.21W$celltype
phenot <- as.character(phenot)
names(phenot) <- rownames(InN.21W@meta.data)
emb <- BGPC.21W @reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './')
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

###pseudotime analysis
library(monocle)
library(data.table)
library(igraph)
library(dplyr)
library(heatmaply)
library(patchwork)
library(gplots)
library(ggplot2)
library(ggsci)
scRNA <- InN_data
scRNA$SampleID <- "SampleID"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("TOF") ] <- "TOF"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
scRNA$SampleID[ as.character(scRNA$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_4") ] <- "Con.20W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_8")] <- "Con.23W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_9") ] <- "Con.25W"
scRNA[['RNA_snn_res.0.5']] <- NULL
head(scRNA@meta.data)
scRNA$celltype_SampleID <- paste0(scRNA$celltype,"*",scRNA$SampleID)
AB <-scRNA[,scRNA@meta.data[["SampleID"]] %in% c("TOF","Con.21W")]
expr_matrix <- as(as.matrix(AB@assays$RNA@counts), 'sparseMatrix')
p_data <- AB@meta.data
p_data$celltype <- AB@active.ident
f_data <- data.frame(gene_short_name = row.names(AB),row.names = row.names(AB))
pd <- new('AnnotatedDataFrame',data = p_data) 
fd <- new('AnnotatedDataFrame',data = f_data) 
cds <- newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,
                      lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes<- row.names(subset(fData(cds), num_cells_expressed >= 10))
diff <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr='~celltype', cores=1)#时间久
head(diff)
deg <- subset(diff, qval < 0.05)
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
write.table(deg,file="InN_21Wmonocle.DEG.xls", col.names=T, row.names=F, sep="\t", quote=F)
ordergene <- row.names(deg)[order(deg$qval)][1:2000]
ordergene
write.csv(ordergene, file = "InN_21W_ordergene.csv", append = FALSE, quote = TRUE, sep = "", row.names = TRUE, col.names = TRUE)
cds.1 <- setOrderingFilter(cds, ordergene)
table(cds.1@featureData@data[["use_for_ordering"]])
pdf("InN_21W_ordergenes.pdf")
plot_ordering_genes(cds.1)
dev.off()
cds.1<- reduceDimension(cds.1, max_components = 2, method = 'DDRTree')
cds.1 <- orderCells(cds.1)
cds.1 <- orderCells(cds.1,root_state = 6)
pdf('InN.21W.monocle+pseudotime.pdf', width = 7, height = 7)
plot_cell_trajectory(cds.1, color_by = 'Pseudotime',size = 1, show_backbone = TRUE)
dev.off()

color2=c(InN_SOX4= "#E5A84B", InN_SGCZ = "#509296",InN_ADARB2 = "#EE7959",InN_RELN = "#9BB496",
         InN_MAF = "#F08080",InN_CCK = "#A64036")
pdf('InN.21W_monocle+celltype.pdf', width = 7, height = 7)
plot_cell_trajectory(cds.1, color_by = 'celltype',size = 1, show_backbone = TRUE)+scale_colour_manual(values = color2)
dev.off()
color2=c(InN_SOX4= "#E5A84B", InN_SGCZ = "#509296",InN_ADARB2 = "#EE7959",InN_RELN = "#9BB496",
         InN_MAF = "#F08080",InN_CCK = "#A64036")
p1 <- plot_complex_cell_trajectory(cds.1, x = 1, y = 2,   color_by = "celltype")+scale_colour_manual(values = color2)
ggsave("InN.21W_monocle+celltype.pdf", plot = p1, width = 10, height = 10)
plot_cell_trajectory(cds.1, color_by = "State",size=1,show_backbone=TRUE)
p1 <- plot_complex_cell_trajectory(cds.1, x = 1, y = 2,   color_by = "State")
ggsave("InN.21W.monocle+State.pdf", plot = p1, width = 10, height = 10)

p1 <- plot_complex_cell_trajectory(cds.1, x = 1, y = 2,   color_by = "SampleID")
ggsave("InN.21W.monocle+SampleID.pdf", plot = p1, width = 10, height = 10)


#SCENIC analysis
library(Seurat)
library(SCENIC)
library(RcisTarget)
library(GENIE3)
library(AUCell)
library(BiocManager)
library(KernSmooth)
library(doParallel)
library(RColorBrewer)
library(SCopeLoomR)
library(data.table)
library(rbokeh)
library(foreach)
library(iterators)
library(parallel)
library(BiocParallel)
library(plumber)
library(Rtsne)
library(stringr)
scenicOptions@fileNames$int
sessionInfo()
SCRNA <-InN_data
SCRNA$SampleID <- "SampleID"
SCRNA$SampleID[ as.character(SCRNA$orig.ident) %in% c("TOF") ] <- "TOF"
SCRNA$SampleID[ as.character(SCRNA$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
SCRNA$SampleID[ as.character(SCRNA$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
SCRNA$SampleID[ as.character(SCRNA$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
SCRNA$SampleID[ as.character(SCRNA$orig.ident) %in% c("Con_4") ] <- "Con.20W"
SCRNA$SampleID[ as.character(SCRNA$orig.ident) %in% c("Con_8")] <- "Con.23W"
SCRNA$SampleID[ as.character(SCRNA$orig.ident) %in% c("Con_9") ] <- "Con.25W"
SCRNA$SampleID  <- factor(SCRNA$SampleID, levels=c("TOF", "Con.21W","Con.22W","Con.17W", "Con.20W", "Con.23W","Con.25W"))
head(SCRNA@meta.data)
SCRNA$celltype_SampleID <- paste0(SCRNA$celltype,"*",SCRNA$SampleID)
AB1 <-SCRNA[,SCRNA@meta.data[["SampleID"]] %in% c("TOF","Con.21W")]
AB1@meta.data[["SampleID"]] <- as.character(AB1@meta.data[["SampleID"]])
table(AB1@meta.data$SampleID)
Idents(AB1) <- "celltype"
save(AB1,file="InN.TOF_21W.RData")
table(AB1@meta.data$celltype)
scenicdata <- subset(AB1,downsample=300)
scenicdata <- RunTSNE(scenicdata,dims.use = 1:20,reduction.use = "pca", dim_embed = 2)
save(scenicdata,file="InN+tsne+scenicdata.RData")
table(scenicdata@active.ident)
exprMat<-GetAssayData(object = scenicdata)#count
dim(exprMat)
exprMat[1:4,1:10]
cellInfo <- scenicdata@meta.data[,c("nCount_RNA","nFeature_RNA","celltype","SampleID","celltype_SampleID")]
head(cellInfo)
cellTypeColumn <- "celltype"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
cbind(table(cellInfo$CellType))
head(cellInfo)
#改颜色
colVars <- list(CellType=c("InN_ADARB2"="#D8B70A",
                           "InN_CCK"="#AD002AB2", 
                           "InN_MAF"="#808180FF", 
                           "InN_RELN"="#E7298A", 
                           "InN_SGCZ"="#9BB496",
                           "InN_SOX4"="#509296" ))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
org="hgnc"
dbDir="/home/lulu/NTD/前端第二次质控/王scenic分析/InN 第二次注释/RcisTarget_databases"
myDatasetTitle="TOF" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,  nCores=10)#初始化
#保存
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
saveRDS(cellInfo, file="int/cellInfo.Rds")
save(exprMat,file = "int/exprMat.Rds")
saveRDS(colVars, file="int/colVars.Rds")
exprMat<-as.matrix(exprMat)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
exprMat_filtered_log <- na.omit(exprMat_filtered_log)
runGenie3(exprMat_filtered_log, scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log,skipTsne = T)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
minPerc <- 1.4
pdf("InN_Scale(1.4).pdf", width = 5, height = 10)
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType_Binarized), center = T, scale=T))
binaryActPerc_subset <- regulonActivity_byCellType_Scaled[which(rowSums(regulonActivity_byCellType_Scaled>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset,
                   color = colorRampPalette(c("#80B9C8", "#FFE3DC", "#ED5A40"))(100),
                   breaks = seq(-1, 2, length.out = 100),
                   border_color = "NA")
dev.off()

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype_SampleID),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

minPerc <- .2
pdf("InN_celltype_SampleID(.2).pdf", width = 10, height = 5)
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
binaryActPerc_subset <- regulonActivity_byCellType_Scaled[which(rowSums(regulonActivity_byCellType>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset,
                   color = colorRampPalette(c("white", "#FFE3DC", "#ED5A40"))(100),
                   breaks = seq(0, 1, length.out = 100),
                   border_color = "NA")
dev.off()
topRegulators <- reshape2::melt(binaryActPerc_subset)
colnames(topRegulators) <- c("Regulon", "celltype_SampleID", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators,"topRegulators_celltypeSampleID.csv")



