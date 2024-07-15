load(obj.combined2.RData)
ExN_data <- obj.combined2[,obj.combined2@meta.data[["celltype"]] %in% c("ExN")]
DefaultAssay(ExN_data) <- "RNA"
ExN_data <-  NormalizeData(ExN_data,normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 2000)
ExN_data <- ScaleData(ExN_data, features = rownames(ExN_data))
ExN_data <- RunPCA(ExN_data, verbose = FALSE)
plot2 <- ElbowPlot(ExN_data,ndims = 50)
dimNums = 10 
ExN_data <- FindNeighbors(ExN_data, reduction = "pca", dims = 1:dimNums)
ExN_data <- FindClusters(ExN_data,resolution = 0.5)
ExN_data <- RunUMAP(ExN_data, dims = 1:dimNums)
DimPlot(ExN_data, reduction = "tsne",cols =c('#B5D465',  '#8EA0C7', '#D790C1', '#EB8F6B', '#82BFA6','#FAD753','#AFBFCF'), pt.size=0.3, label = TRUE,repel = TRUE)
ExN_makers=c('TUBB3', "FOXP2",	"SLC17A7",	"ADRA2A",	"CARTPT",	"CNR1",	"CTGF",	"CUX2",	"DCX",	"EPHA6",	"ETV1",	"FAM84B",	"FEZF2",	"GLRA3",	"GRIK4",	"HTR2C",	"LAMP5",	"MEIS2",	"NNAT",	"NR4A2",	"NTING2",	"NXPH4",	"OPRK1",	"OTOF",	"PCDH8",	"PCP4",	"PRSS12",	"RASGRF2",	"RBFOX1",	"RORB",	"RPRM",	"RSPO1",	"RXFP1",	"SLA2",	"SULF1",	"SYT6",	"THEMIS",	"TLE4",	"TOX",	"TPBG",	"UNC5D",	"XIST")
DotPlot(ExN_data, assay = "RNA", features = ExN_makers, cols = c("lightgrey", "red"), col.min = 0,col.max = 2) + theme(axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 15))
ExN_data <- ExN_data[,ExN_data@meta.data[["seurat_clusters"]] %in% c('0',"1",	"2",	"3", "5")]
DefaultAssay(ExN_data) <- "RNA" 
ExN_data <-  NormalizeData(ExN_data,normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 2000)
ExN_data <- ScaleData(ExN_data, features = rownames(ExN_data))
ExN_data <- RunPCA(ExN_data, verbose = FALSE)
plot2 <- ElbowPlot(ExN_data,ndims = 50)
dimNums = 10 
ExN_data <- FindNeighbors(ExN_data, reduction = "pca", dims = 1:dimNums)
ExN_data <- FindClusters(ExN_data,resolution = 0.5)
ExN_data <- RunUMAP(ExN_data, dims = 1:dimNums)
p9 <- DimPlot(ExN_data, cols = c('#FFFFC4',  '#D8E7D2', '#B1635D','#98C0E4', '#8EC9CF', '#9C96BE','#E7AEA4'), reduction = "umap", group.by = "ident", pt.size=0.1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "ExN_umap_tsne1_1.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')
Makers <- FindAllMarkers(ExN_data, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)


##subtype annotation
current.cluster.ids <- c("0",	"1",	"2", "3",	"4",	"5",	"6")
new.cluster.ids <- c("ExN_TLE4", "ExN_CNTNAP2",	"ExN_ENC1",	"ExN_CNTN5",	"ExN_TUBB3",	"ExN_MEIS2", "ExN_POU6F2") 
ExN_data@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(ExN_data@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(ExN_data@meta.data$celltype)
current.cluster.ids <- c("0",	"1",	"2", "3",	"4",	"5",	"6")
names(new.cluster.ids) <- levels(ExN_data)
ExN_data <- RenameIdents(ExN_data, new.cluster.ids)
DimPlot(ExN_data, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
save(ExN_data,file="ExN_data+subtype.RData")

llcolour=c('#FFFFC4',  '#D8E7D2', '#B1635D','#98C0E4', '#8EC9CF', '#9C96BE','#E7AEA4')
umap=ExN_data@reductions$umap@cell.embeddings%>%as.data.frame()%>%cbind(celltype=ExN_data@meta.data$celltype)
p9 <-ggplot(umap,aes(x=UMAP_1,y=UMAP_2,fill=celltype))+geom_point(size=3.5,colour='grey35',shape=21)+scale_fill_manual(values=allcolour)+theme_classic()
ggsave(filename = "ExN+subtype+UMAP.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')


library(pheatmap)
DefaultAssay(ExN_data) <- "RNA"
markers <- FindAllMarkers(ExN_data, logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
head(markers)
sig_markers <- markers %>% group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
head(sig_markers)
genes<- unique(sig_markers$gene)
aver_dt <- AverageExpression(ExN_data,features = genes,group.by ='celltype',slot="data")
aver_dt<-as.data.frame(aver_dt$RNA)
aver_dt[1:6,1:6]
cell_anno<-data.frame(cell_anno=colnames(aver_dt),row.names =colnames(aver_dt))
celltype_col<-c('#FFFFC4',  '#D8E7D2', '#B1635D','#98C0E4', '#8EC9CF', '#9C96BE','#E7AEA4')
names(celltype_col)<-cell_anno$cell_anno
anno_col<-list(cell_anno=celltype_col)
anno_col
p1=pheatmap(as.matrix(aver_dt),scale='row',cellwidth = 10, cellheight = 5,cluster_rows = FALSE,cluster_cols =FALSE,annotation_col = cell_anno,annotation_colors = anno_col,border_color = 'white')
ggsave(filename = "ExN_data+TOP10+heatmap.pdf", plot = p1, device = 'pdf', width = 10, height = 18, units = 'cm')

ExN_data$SampleID <- "SampleID"
ExN_data$SampleID[ as.character(ExN_data$orig.ident) %in% c("TOF") ] <- "TOF"
ExN_data$SampleID[ as.character(ExN_data$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
ExN_data$SampleID[ as.character(ExN_data$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
ExN_data$SampleID[ as.character(ExN_data$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
ExN_data$SampleID[ as.character(ExN_data$orig.ident) %in% c("Con_4") ] <- "Con.20W"
ExN_data$SampleID[ as.character(ExN_data$orig.ident) %in% c("Con_8")] <- "Con.23W"
ExN_data$SampleID[ as.character(ExN_data$orig.ident) %in% c("Con_9") ] <- "Con.25W"
AB <- ExN_data[,ExN_data@meta.data[["SampleID"]] %in% c("TOF","Con.21W")]
p9 <- DimPlot(AB, cols = c('#D9E6B7','#8EA0C7'), reduction = "umap", group.by = "SampleID", pt.size=0.5, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "ExNUMAP+TOF_Con21.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')
p9 <- DimPlot(ExN_data, cols = c("#8DD3C7", "#FDB462", "#BEBADA", "#FB8072", "#80B1D3","#FFFFB3", "#B3DE69"), reduction = "umap", group.by = "SampleID", pt.size=0.01, label = F,label.size=5, repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "ExN+umap+SampleID.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')
AB <- ExN_data[,ExN_data@meta.data[["SampleID"]] %in% c("TOF")]
p9 <- DimPlot(AB, cols = c('#FFFFC4',  '#D8E7D2', '#B1635D','#98C0E4', '#8EC9CF', '#9C96BE','#E7AEA4'), reduction = "umap", group.by = "celltype", pt.size=0.5, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "ExNUMAP+TOF.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')

AB <- ExN_data[,ExN_data@meta.data[["SampleID"]] %in% c("Con.21W")]
p10 <- DimPlot(AB, cols = c('#FFFFC4',  '#D8E7D2', '#B1635D','#98C0E4', '#8EC9CF', '#9C96BE','#E7AEA4'), reduction = "umap", group.by = "celltype", pt.size=0.5, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "ExNUMAP+Con.21W.pdf", plot = p10, device = 'pdf', width = 15, height = 15, units = 'cm')


###finding DEGs
table(ExN_data@meta.data$orig.ident)
ExN_data$group=str_replace(ExN_data$orig.ident,"_//d","")
table(ExN_data@meta.data$group)
table(ExN_data@meta.data$celltype)
Idents(ExN_data)='group'
levels(ExN_data) 
group2 = c('TOF',"Con_21","Con_21","Con_22","Con_20", "Con_22", 'Con_17', 'Con_17', 'Con_23','Con_25')
names(group2) = levels(ExN_data)
ExN_data <- RenameIdents(ExN_data, group2)
levels(ExN_data)
ExN_data$group2 <- Idents(ExN_data)
ExN_data$celltype.group2 <- paste(ExN_data$celltype,ExN_data$group2, sep = "_") ##细胞类型前面加上分组信息
Idents(ExN_data) <- "celltype.group2"
table(ExN_data@active.ident)

ExN.TOF_CON_21.diff <- FindMarkers(ExN_data, ident.1 = "ExN_TUBB3_TOF", ident.2 = "ExN_TUBB3_Con_21", logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
write.table(ExN.TOF_CON_21.diff,file="21W_ExN.TUBB3.TOF_CON.diff_genes.txt",quote=F,sep="\t",row.names=T,col.names=T)

ExN.TOF_CON.diff <- FindMarkers(ExN_data, ident.1 = "ExN_CNTNAP2_TOF", ident.2 = "ExN_CNTNAP2_Con_21", logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
write.table(ExN.TOF_CON.diff,file="_21W_ExN.CNTNAP2.TOF_CON.diff_genes.txt",quote=F,sep="\t",row.names=T,col.names=T)

ExN.TOF_CON_21.diff <- FindMarkers(ExN_data, ident.1 = "ExN_TLE4_TOF", ident.2 = "ExN_TLE4_Con_21", logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
write.table(ExN.TOF_CON_21.diff,file="21W_ExN.TLE4.TOF_CON.diff_genes.txt",quote=F,sep="\t",row.names=T,col.names=T)

ExN.TOF_CON_21.diff <- FindMarkers(ExN_data, ident.1 = "ExN_CNTN5_TOF", ident.2 = "ExN_CNTN5_Con_21", logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
write.table(ExN.TOF_CON_21.diff,file="21W_ExN.CNTN5.TOF_CON.diff_genes.txt",quote=F,sep="\t",row.names=T,col.names=T)

ExN.TOF_CON_21.diff <- FindMarkers(ExN_data, ident.1 = "ExN_MEIS2_TOF", ident.2 = "ExN_MEIS2_Con_21", logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
write.table(ExN.TOF_CON_21.diff,file="21W_ExN.MEIS2.TOF_CON.diff_genes.txt",quote=F,sep="\t",row.names=T,col.names=T)


ExN.TOF_CON_21.diff <- FindMarkers(ExN_data, ident.1 = "ExN_POU6F2_TOF", ident.2 = "ExN_POU6F2_Con_21", logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
write.table(ExN.TOF_CON_21.diff,file="21W_ExN.POU6F2.TOF_CON.diff_genes.txt",quote=F,sep="\t",row.names=T,col.names=T)


ExN.TOF_CON_21.diff <- FindMarkers(ExN_data, ident.1 = "ExN_ENC1_TOF", ident.2 = "ExN_ENC1_Con_21", logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
write.table(ExN.TOF_CON_21.diff,file="21W_ExN.ENC1.TOF_CON.diff_genes.txt",quote=F,sep="\t",row.names=T,col.names=T)


#Mfuzz analysis
library(Mfuzz)
library(CellChat)
library(Seurat)
library(progeny)
library(tidyr)
library(tibble)
library(dplyr)
library(viridisLite)
library(ggplot2)
invisible(utils::memory.limit(100000))
library("BiocParallel")
register(SnowParam(6))
Sys.setenv("VROOM_CONNECTION_SIZE"=88888888)
scRNA<- ExN_data
scRNA$SampleID <- "SampleID"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("TOF") ] <- "TOF"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
scRNA$SampleID[ as.character(scRNA$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_4") ] <- "Con.20W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_8")] <- "Con.23W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_9") ] <- "Con.25W"
scRNA1 <- scRNA[,scRNA@meta.data[["SampleID"]] %in% c("Con.17W","Con.20W", "Con.21W","Con.22W","Con.23W","Con.25W")]
scRNA1=subset(scRNA1,features= rownames(scRNA1@assays$RNA@scale.data))
scRNA1@meta.data[["SampleID"]]<-factor(scRNA1@meta.data[["SampleID"]], levels=c("Con.17W","Con.20W", "Con.21W","Con.22W","Con.23W","Con.25W"))
scRNA1=ScaleData(scRNA1)
age.averages <- AverageExpression(scRNA1,group.by = "SampleID")
df1 <- age.averages[["RNA"]]
mat <- as.matrix(df1)
dt <- new("ExpressionSet",exprs = mat)
dim(dt)
dt.r <- filter.NA(dt, thres=0.25)
dt.f <- fill.NA(dt.r,mode="mean")
tmp <- filter.std(dt.f,min.std=0.05)
dt.s <- standardise(tmp)
df.s <- dt.s@assayData$exprs
m1 <- mestimate(dt.s)
set.seed(007)
cl <- mfuzz(dt.s,c=6,m=m1)
mfuzz.plot(dt.s,cl,mfrow=c(2,3), new.window= FALSE, time.labels=colnames(dt.s))
library(RColorBrewer)
mycol <- c("#E5A84B","#EE7959","#8EC9CF")
mycolor <- colorRampPalette(mycol)(20)
pdf("RNA_Con_ExN expression changes trends_0.05.pdf", width = 12,height = 10)
mfuzz.plot(dt.s,cl,mfrow=c(3,3), new.window= FALSE, time.labels=colnames(dt.s), colo = mycolor)
dev.off()
gene_cluster <- data.frame(df.s,cluster=cl$cluster)
head(gene_cluster)
write.csv(gene_cluster, 'RNA_Con_ExN expression changes trends_0.05.csv')
#cluster1
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix)
library(readxl)
library(ggpubr)


retina_gene<-read_xlsx("cluster 1.xlsx")
gene<-as.list(retina_gene)
AB<-AddModuleScore(scRNA1, features = gene, ctrl = 100, name = "retina")
colnames(AB@meta.data)[8]<-"retina_Score"
AB@meta.data[["celltype"]]<-factor(AB@meta.data[["celltype"]], levels=c("ExN_MEIS2","ExN_CNTN5",	"ExN_CNTNAP2",	"ExN_POU6F2",	"ExN_TLE4", "ExN_TUBB3","ExN_ENC1"))
library(ggpubr)
color2=c(ExN_CNTN5='#FFFFC4', ExN_CNTNAP2='#D8E7D2', ExN_MEIS2='#B1635D',ExN_POU6F2='#98C0E4', ExN_SOX11='#8EC9CF', ExN_TLE4='#9C96BE',ExN_TUBB3='#E7AEA4')
p.AddModuleScore <- ggviolin(AB@meta.data, x = "celltype", y = "retina_Score",
                             color = "celltype",add = 'mean_sd',fill = 'celltype',
                             add.params = list(color = "black")) + 
  scale_color_manual(values = color2) + 
  scale_fill_manual(values = color2) +
  NoLegend() + labs(x = '')
p.AddModuleScore
ggsave('RNA_Con_ExN expression changes trends_0.05_cluster1.pdf',width = 10,height = 5)

##INTERSECTING GENE
library(ggplot2)
library(ggsci)
library(sf)
library(ggVennDiagram)
CNTNAP2 <- read_excel("CNTNAP2.xlsx")
MEIS2 <- read_excel("MEIS2.xlsx")
POU6F2 <- read_excel("POU6F2.xlsx")
TLE4<- read_excel("TLE4.xlsx")
gene.list = list()
gene.list <- list(CNTN5 = CNTN5$GENE,
                  CNTNAP2=CNTNAP2$GENE,
                  MEIS2=MEIS2$GENE,
                  POU6F2=POU6F2$GENE,
                  TLE4=TLE4$GENE)
pdf("ExN5subtype+shareddown+venny.pdf",width = 8,height = 7)
ggVennDiagram(gene.list[1:5], set_size = 8,label_alpha=0,label_size =6,
              edge_size = 1.2,label ="count") +
  scale_color_lancet()+
  scale_fill_gradient(low="gray100",high = color1,guide="none")
dev.off()

data <- read.csv('shared down in 5subtype.txt',sep= '')
data1 <- read_excel("cluster 1.xlsx")
gene.list = list()
gene.list <- list(downin5 = data$GENE,
                  cluster1=data1$GENE)
summary(gene.list) 
pdf("venny(downin5andcluster1).pdf",width = 8,height = 7)
ggVennDiagram(gene.list[1:2], set_size = 8,label_alpha=0,label_size =6,
              edge_size = 1.2,label ="count") +
  scale_color_lancet()+
  scale_fill_gradient(low="gray100",high = color2,guide="none")
dev.off()

Data <- read.csv('37intersectinggene.txt',sep= '')
mygene <- Data$gene
gene.df <- bitr(mygene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)

ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoBP@result,'37intersectinggeneGO.BP.result.csv')



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
scRNA<- ExN_data
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
rm(scRNA)
DefaultAssay(AB)<- "RNA"
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
ordergene <- row.names(deg)[order(deg$qval)][1:2000]
ordergene
write.csv(ordergene, file = "ordergene.csv", append = FALSE, quote = TRUE, sep = "", row.names = TRUE, col.names = TRUE)
cds <- setOrderingFilter(cds, ordergene)
table(cds@featureData@data[["use_for_ordering"]])
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = 1)
Time_diff <-differentialGeneTest(cds[ordergene,], cores = 1,
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff,"Time_diff_all.csv", row.names = F)
topNGenes <- top_n(Time_diff, n = 500, desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()

pht <- plot_pseudotime_heatmap(
  cds[topNGenes,],
  num_clusters = 5,
  show_rownames = T,
  return_heatmap = T
)
pht
ggsave(pht, file = "time_diff_top500.pdf", width = 5, height = 20)
cds.1 <- orderCells(cds,root_state = 4)
colour = c("#3A5A7D", "#82B289")
pdf('monocle+pseudotime.pdf', width = 7, height = 7)
plot_cell_trajectory(cds.1, color_by = 'Pseudotime',size = 1, show_backbone = TRUE)
dev.off()
pdf('monocle+celltype.pdf', width = 7, height = 7)
plot_cell_trajectory(cds.1, color_by = 'celltype',size = 1, show_backbone = TRUE)
dev.off()
pdf("monocle+state.pdf",width = 7,height = 7)
plot_cell_trajectory(cds.1, color_by = "State",size=1,show_backbone=TRUE)
dev.off()
p7 <- plot_cell_trajectory(cds.1, color_by = "SampleID",size=1,show_backbone=TRUE)
ggsave(p7, file = "monocle+SampleID.pdf",width = 10, height = 10)
pdf('monocle+celltype_facet.pdf', width = 7, height = 14)
plot_cell_trajectory(cds.1, color_by = 'celltype',size = 1, show_backbone = TRUE)+ facet_wrap("~celltype", nrow = 4)
dev.off()
color <-c( '#FFFFC4',  '#D8E7D2', '#B1635D','#98C0E4', '#8EC9CF', '#9C96BE','#E7AEA4')
pdf("monocle+celltype.pdf",width = 7,height = 7)
plot_cell_trajectory(cds.1, color_by = "celltype",size=1,show_backbone=TRUE)+scale_colour_manual(values = color)
dev.off()
color <-c( "#82B289","#82B289",'#D8E7D2', '#D8E7D2','#B1635D','#B1635D',
           '#8EC9CF','#8EC9CF', "#16ABD0","#16ABD0",'#FFFFC4', '#FFFFC4', '#9C96BE','#9C96BE')
pdf('monocle+celltype_SampleID_facet.pdf', width = 10, height = 35)
plot_cell_trajectory(cds.1, color_by = 'celltype_SampleID',size = 1, show_backbone = TRUE)+ facet_wrap("~celltype_SampleID", nrow = 7)+scale_colour_manual(values = color)
dev.off()

Time_diff <-differentialGeneTest(cds.1[ordergene,], cores = 1,
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
pdf("Time_heatmap_select.pdf", width = 5, height = 3)
marker_genes <- row.names(subset(fData(cds.1),
                                 gene_short_name %in% c("FABP7","STMN1","NEUROD6","STMN2","FLRT2","SOX11","TLE4","CNTN5","POU6F2","TUBB3","CNTNAP2","MEIS2")))
diff_test_res <- differentialGeneTest(cds.1[marker_genes,],fullModelFormulaStr ="~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
plot_pseudotime_heatmap(cds.1[sig_gene_names,], num_clusters = 1, cores = 1, show_rownames = T,hmcols = colorRampPalette(c("#83AED0","#E94A19"))(500))
dev.off()

BEAM_res <- BEAM(cds.1[ordergene,], branch_point = 2, cores = 1,progenitor_method = "duplicate") 

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res, "BEAM_res.csv", row.names = F)
BEAM_gene_selected <- row.names(subset(fData(cds.1),
                                       gene_short_name %in% c("SOX11","EPHA3","ENC1","CNR1","POU3F2",
                                                              "POU6F2","SYT4","KIT")))
p <- plot_genes_branched_heatmap(cds.1[BEAM_gene_selected,],  branch_point = 2,cluster_rows = F,
                                 num_clusters = 3,show_rownames = T, return_heatmap = T, hmcols = colorRampPalette(c("#83AED0","#E94A19"))(60))
ggsave("ExN.21W.branch_heatmaptop_selected1.pdf", p$ph_res, width = 7, height = 7)


devtools::install_github("junjunlab/ClusterGVis")
library(ClusterGVis)
BEAM_gene_selected <- row.names(subset(fData(cds.1),
                                       gene_short_name %in% c("SOX11","EPHA3","ENC1","CNR1","POU3F2",
                                                              "SYT4","KIT")))
markGenes <- c("SOX11","EPHA3","ENC1","CNR1","POU3F2",
               "SYT4","KIT")
df <- plot_genes_branched_heatmap2(cds.1[BEAM_gene_selected,],
                                   branch_point = 2,
                                   num_clusters = 1,
                                   cores = 1,
                                   use_gene_short_name = T,
                                   show_rownames = T)
str(df)
visCluster(object = df,plot.type = "heatmap")

visCluster(object = df,plot.type = "heatmap",
           pseudotime_col = c("#82B289","grey","#3A5A7D"))
pdf(file = "heatmap.pdf",height = 6,width = 8)
visCluster(object = df,plot.type = "both",  line.side = "right",  markGenes = markGenes,
           markGenes.side = "left", pseudotime_col = c("#82B289","grey","#3A5A7D"),
           ht.col.list =list(col_range = c(-2, 0, 2),col_color = c("#16ABD0", "#FFE3DC", "#ED5A40")),
           ctAnno.col = c("#C389BC", "#C98F8D", "#78ABDB","#82B289"))
dev.off()


keygenes <- c("FLRT2")
cds_subset <- cds.1[keygenes,]
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "celltype")
p4 <- plot_genes_in_pseudotime(cds_subset, color_by = "SampleID")
plotc <- p1|p2|p3|p4
ggsave("ExN.21W.keygenes.FLRT2.pdf", plot = plotc, width = 20, height = 4)

pdf("ExN.21W.keygene.pdf",width = 4,height = 4)
genes <- row.names(subset(fData(cds.1),
                          gene_short_name %in% c("EPHA3","SYT4","CNR1")))
plot_genes_branched_pseudotime(cds.1[genes,],
                               branch_point = 2,
                               color_by = "SampleID",
                               ncol = 1)
dev.off()


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
library(ggplot2)
library(Rtsne)
scenicOptions@fileNames$int
sessionInfo()
install.packages("VlnPlot")

scRNA <-ExN_data
###合并某些接近的周数样本
###样本分组
scRNA$SampleID <- "SampleID"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("TOF") ] <- "TOF"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
scRNA$SampleID[ as.character(scRNA$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_4") ] <- "Con.20W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_8")] <- "Con.23W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_9") ] <- "Con.25W"
scRNA$SampleID  <- factor(scRNA$SampleID, levels=c("TOF", "Con.21W","Con.22W","Con.17W", "Con.20W", "Con.23W","Con.25W"))
head(scRNA@meta.data)
scRNA$celltype_SampleID <- paste0(scRNA$celltype,"*",scRNA$SampleID)
AB <-scRNA[,scRNA@meta.data[["SampleID"]] %in% c("TOF","Con.21W")]
AB@meta.data[["SampleID"]] <- as.character(AB@meta.data[["SampleID"]])
AB@meta.data[["celltype"]] <- as.character(AB@meta.data[["celltype"]])
table(AB@meta.data$SampleID)
Idents(AB) <- "celltype"
save(AB,file="ExN.TOF_21W.RData")
table(AB@meta.data$celltype)
table(scenicdata@meta.data$SampleID)

scenicdata <- subset(AB,downsample=300)
scenicdata <- RunTSNE(scenicdata,dims.use = 1:20,reduction.use = "pca", dim_embed = 2)
save(scenicdata,file="ExN+scenicdata.RData")

exprMat<-GetAssayData(object = scenicdata)#count
dim(exprMat)
exprMat[1:4,1:10]
cellInfo <- scenicdata@meta.data[,c("nCount_RNA","nFeature_RNA","celltype","SampleID","celltype_SampleID")]#创建cellInfo
head(cellInfo)
cellTypeColumn <- "celltype"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
cbind(table(cellInfo$CellType))
head(cellInfo)
colVars <- list(CellType=c("ExN_CNTN5"="skyblue",
                           "ExN_CNTNAP2"="darkblue",
                           "ExN_MEIS2"="peru",
                           "ExN_POU6F2"="salmon", 
                           "ExN_ENC1"="#E31A1C", 
                           "ExN_TLE4"="#FDBF6F",
                           "ExN_TUBB3"="#FF7F00"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
org="hgnc"
dbDir="/home/lulu/NTD/前端第二次质控/王scenic分析/ExN/cisTarget_databases"
myDatasetTitle="TOF" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir,  nCores=10)
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
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype_SampleID),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
minPerc <- .3
pdf("ExN_celltype_SampleID(0.3).pdf", width = 9, height = 4)
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
binaryActPerc_subset <- regulonActivity_byCellType_Scaled[which(rowSums(regulonActivity_byCellType>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset,
                   color = colorRampPalette(c("white","pink","red"))(100),
                   breaks = seq(-1, 2, length.out = 100),
                   border_color = "NA")
dev.off()
