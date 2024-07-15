load("obj.combined2.RData")
GC_data <- obj.combined2[,obj.combined2@meta.data[["celltype"]] %in% c("GC")]
DefaultAssay(GC_data) <- "RNA" 
colnames(GC_data@meta.data)
GC_data[['seurat_clusters']] <- NULL
GC_data[['integrated_snn_res.0.6']] <- NULL
GC_data[['celltype']] <- NULL
GC_data <-  NormalizeData(GC_data,normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 2000)
GC_data <- ScaleData(GC_data, features = rownames(GC_data))
GC_data <- RunPCA(GC_data, verbose = FALSE)
plot2 <- ElbowPlot(GC_data,ndims = 50)
dimNums = 10
GC_data <- FindNeighbors(GC_data, reduction = "pca", dims = 1:dimNums)
C_data <- FindClusters(GC_data,resolution = 0.6)
GC_data <- RunUMAP(GC_data, dims = 1:dimNums)
DefaultAssay(GC_data) <- "integrated"
markers <- FindAllMarkers(GC_data, logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
GC_makers=c('TNC','CSPG5','SOX2',"ANGPT1","EGFR","ASCL1","PDGFRA","OLG1","OLG2",'CSPG4',	"PCDH15",	"PMP2",	"COL20A1","GRID2",	"LSAMP", 'INPP5D','CSF1R', "SPI",	"PTRRC", 'GLUL',	"GJA1",	'SOX9', 'AQP4', 'SLC1A3','SLC1A2', "NDRG2",'SPP1','MAF','C3', 'CSF3R','CX3CR1', 'FLT1','CLDN5','VWF','PDGFRB','ABCC9')
pdf('GC_maker_exp_dotplot.pdf',width =15,height = 5)
DotPlot(GC_data, assay = "RNA", features = GC_makers, cols = c("lightgrey", "red"), col.min = 0,col.max = 2) + theme(axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size = 15))
dev.off()
current.cluster.ids <- c("0",	"1",	"2", "3",	"4",	"5",	"6",'7','8')
new.cluster.ids <- c("BGPC","BGPC",	"BGPC",	"BGPC",	"Endo",	"Endo", "Peri",'Endo','MG') 
GC_data@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(GC_data@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(GC_data@meta.data$celltype)
current.cluster.ids <- c("0",	"1",	"2", "3",	"4",	"5",	"6",'7','8')
names(new.cluster.ids) <- levels(GC_data)
GC_data <- RenameIdents(GC_data, new.cluster.ids)
DimPlot(GC_data, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
p9 <- DimPlot(GC_data, cols = c('#B5D465',  '#8EA0C7', '#EB8F6B', '#D790C1'), reduction = "umap", group.by = "ident", pt.size=1, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "UMAP.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')
GC_data$SampleID <- "SampleID"
GC_data$SampleID[ as.character(GC_data$orig.ident) %in% c("TOF") ] <- "TOF"
GC_data$SampleID[ as.character(GC_data$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
GC_data$SampleID[ as.character(GC_data$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
GC_data$SampleID[ as.character(GC_data$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
GC_data$SampleID[ as.character(GC_data$orig.ident) %in% c("Con_4") ] <- "Con.20W"
GC_data$SampleID[ as.character(GC_data$orig.ident) %in% c("Con_8")] <- "Con.23W"
GC_data$SampleID[ as.character(GC_data$orig.ident) %in% c("Con_9") ] <- "Con.25W"
table(GC_data@meta.data$SampleID)
AB <- GC_data[,GC_data@meta.data[["SampleID"]] %in% c("TOF","Con.21W")]
p9 <- DimPlot(AB, cols = c('#8EA0C7','#EB8F6B'), reduction = "umap", group.by = "SampleID", pt.size=2.5, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "GC+TOF_Con21.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')
AB <- GC_data[,GC_data@meta.data[["SampleID"]] %in% c("TOF")]
p9 <- DimPlot(AB, cols = c('#B5D465',  '#8EA0C7', '#EB8F6B', '#D790C1'), reduction = "umap", group.by = "celltype", pt.size=2.5, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "GC+TOF.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')
AB <- GC_data[,GC_data@meta.data[["SampleID"]] %in% c("Con.21W")]
p10 <- DimPlot(AB, cols = c('#B5D465',  '#8EA0C7', '#EB8F6B', '#D790C1'), reduction = "umap", group.by = "celltype", pt.size=2.5, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "GC+Con.21W.pdf", plot = p10, device = 'pdf', width = 15, height = 15, units = 'cm')
p9 <- DimPlot(GC_data, cols = c("#8DD3C7", "#FDB462", "#BEBADA", "#FB8072", "#80B1D3","#FFFFB3", "#B3DE69"), reduction = "umap", group.by = "SampleID",pt.size=2.5, label = F,label.size=5, repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "GC+umap+SampleID.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')

DefaultAssay(GC_data) <- "RNA"
markers <- FindAllMarkers(GC_data, logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
write.table(markers,file="GC_all_diff_marker_genes.txt",quote=F,sep="\t",row.names=F,col.names=T)

#BGPC+GOBP analysis
Data <- read.csv('BGPC.txt',sep= '')
mygene <- Data$gene 
gene.df <- bitr(mygene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)

ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoBP@result,'BGPC+GO.BP.result.csv')


#BGPC subtype analysis
BGPC <- GC_data[,GC_data@meta.data[["celltype"]] %in% c("BGPC")]
DefaultAssay(BGPC) <- "RNA" 
colnames(BGPC@meta.data)
BGPC[['seurat_clusters']] <- NULL
BGPC[['RNA_snn_res.0.6']] <- NULL
BGPC[['celltype']] <- NULL

BGPC <-  NormalizeData(BGPC,normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 2000)
BGPC <- ScaleData(BGPC, features = rownames(BGPC))
BGPC <- RunPCA(BGPC, verbose = FALSE)
plot2 <- ElbowPlot(BGPC,ndims = 50)
dimNums = 8 
BGPC <- FindNeighbors(BGPC, reduction = "pca", dims = 1:dimNums)
BGPC <- FindClusters(BGPC,resolution = 0.6)
BGPC <- RunUMAP(BGPC, dims = 1:dimNums)
p9 <- DimPlot(BGPC, cols = c('#B5D465',  '#8EA0C7', '#D790C1','#EB8F6B', '#82BFA6','#FAD753','#AFBFCF','#845EC2', '#D65DB1'), reduction = "umap", group.by = "ident", pt.size=2, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
Makers <- FindAllMarkers(BGPC, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.25)
current.cluster.ids <- c("0",	"1",	"2", "3",	"4",	"5")
new.cluster.ids <- c("GPC",	"high_NRG1",	"GPC", "OL",	"OL",	"high_MKI67") 
BGPC@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(BGPC@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(BGPC@meta.data$celltype)
current.cluster.ids <- c("0",	"1",	"2", "3",	"4",	"5")
names(new.cluster.ids) <- levels(BGPC)
BGPC<- RenameIdents(BGPC, new.cluster.ids)
p9 <- DimPlot(BGPC, cols = c('#B5D465',  '#8EA0C7','#D790C1','#EB8F6B'), reduction = "umap", group.by = "celltype", pt.size=2, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "BGPC+umap.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')
BGPC$SampleID <- "SampleID"
BGPC$SampleID[ as.character(BGPC$orig.ident) %in% c("TOF") ] <- "TOF"
BGPC$SampleID[ as.character(BGPC$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
BGPC$SampleID[ as.character(BGPC$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
BGPC$SampleID[ as.character(BGPC$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
BGPC$SampleID[ as.character(BGPC$orig.ident) %in% c("Con_4") ] <- "Con.20W"
BGPC$SampleID[ as.character(BGPC$orig.ident) %in% c("Con_8")] <- "Con.23W"
BGPC$SampleID[ as.character(BGPC$orig.ident) %in% c("Con_9") ] <- "Con.25W"

p9 <- DimPlot(BGPC, cols = c("#8DD3C7", "#FDB462", "#BEBADA", "#FB8072", "#80B1D3","#FFFFB3", "#B3DE69"), reduction = "umap", group.by = "SampleID",pt.size=2, label = F,label.size=5, repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "BGPC+umap+SampleID.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')

AB <- BGPC[,BGPC@meta.data[["SampleID"]] %in% c("TOF","Con.21W")]
p9 <- DimPlot(AB, cols = c('#8EA0C7','#EB8F6B'), reduction = "umap", group.by = "SampleID", pt.size=2.5, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "BGPC+TOF_Con21.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')

AB <- BGPC[,BGPC@meta.data[["SampleID"]] %in% c("TOF")]
p9 <- DimPlot(AB, cols = c('#B5D465',  '#8EA0C7','#D790C1','#EB8F6B'), reduction = "umap", group.by = "celltype", pt.size=2, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "BGPC+TOF.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')
AB <- BGPC[,BGPC@meta.data[["SampleID"]] %in% c("Con.21W")]
p9 <- DimPlot(AB, cols = c('#B5D465',  '#8EA0C7','#D790C1','#EB8F6B'), reduction = "umap", group.by = "celltype", pt.size=2, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "BGPC+Con.21W.pdf", plot = p9, device = 'pdf', width = 15, height = 15, units = 'cm')

DefaultAssay(BGPC) <- "RNA"
celltype_umap<- FeaturePlot(BGPC, reduction = "umap", features = c('SPARCL1','NRG1','PCDH15','MKI67'),ncol = 2,pt.size = 2,cols =c("lightgrey", "#CD0000")) 
ggsave(filename =  "BGPC+marker_umap.pdf", plot = celltype_umap, device = 'pdf', width = 40, height = 40, units = 'cm')
genes_to_check = c("APOE","VIM",'SPARCL1',"NRG1","MKI67","TOP2A","PCDH15","PDGFRA")
library(stringr)  
genes_to_check=str_to_upper(unique(genes_to_check))
genes_to_check

th=theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
p_all_markers <- DotPlot(BGPC, features = genes_to_check,
                         assay='RNA' ,group.by = 'celltype' )  + coord_flip()+th
p_all_markers

data<-p_all_markers$data

colnames(data)

colnames(data)<-c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")

unique(data$`Precent Expressed`)

p1 = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "#498EA4", high = "#E54924")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45, 
                                    vjust = 0.5, hjust=0.5),
        axis.text.y  = element_text(color="black",size=12),
        legend.text = element_text(size =12,color="black"),
        legend.title = element_text(size =12,color="black"),
        axis.title.y=element_text(vjust=1,  
                                  size=16)
  )+labs(x=" ",y = "Features")
p1
ggsave(file="maker_exp_dotplot.pdf",plot = p1,height= 6,width = 5)

#GENE EXPRESSION OF ERBB3 and PDGFRA 
AB <- BGPC[,BGPC@meta.data[["SampleID"]] %in% c("TOF")]
celltype_umap<- FeaturePlot(AB, reduction = "umap", features = c('ERBB3'),ncol = 1,pt.size = 2,cols =c("lightgrey", "#CD0000")) 
ggsave(filename =  "ERBB3+TOF+marker_umap.pdf", plot = celltype_umap, device = 'pdf', width = 15, height = 15, units = 'cm')
celltype_umap<- FeaturePlot(AB, reduction = "umap", features = c('PDGFRA'),ncol = 1,pt.size = 2,cols =c("lightgrey", "#CD0000")) 
ggsave(filename =  "PDGFRA+TOF+marker_umap.pdf", plot = celltype_umap, device = 'pdf', width = 15, height = 15, units = 'cm')
AB <- BGPC[,BGPC@meta.data[["SampleID"]] %in% c("Con.21W")]
celltype_umap<- FeaturePlot(AB, reduction = "umap", features = c('ERBB3'),ncol = 1,pt.size = 2,cols =c("lightgrey", "#CD0000")) 
ggsave(filename =  "ERBB3+Con.21W+marker_umap.pdf", plot = celltype_umap, device = 'pdf', width = 15, height = 15, units = 'cm')
celltype_umap<- FeaturePlot(AB, reduction = "umap", features = c('PDGFRA'),ncol = 1,pt.size = 2,cols =c("lightgrey", "#CD0000")) 
ggsave(filename =  "PDGFRA+Con.21W+marker_umap.pdf", plot = celltype_umap, device = 'pdf', width = 15, height = 15, units = 'cm')

#TOP20 HEATMAP
library(ClusterGVis)
library(org.Hs.eg.db)
library(Seurat)
library(SeuratObject)
markers.all <- Seurat::FindAllMarkers(BGPC,
                                      only.pos = TRUE,
                                      min.pct = 0.15,
                                      logfc.threshold = 0.25)
markers <- markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)
head(markers)
markers1 <- markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)
st.data <- prepareDataFromscRNA(object = BGPC,
                                diffData = markers1,
                                showAverage = TRUE)
str(st.data)
markGenes = c("SPARCL1","ID2","PCDH15","PDGFRA","NRG1","MKI67","TOP2A")
# line plot
P2<- visCluster(object = st.data,
                plot.type = "line",mline.col="#C34A36",mline.size=1,ncol=1)
ggsave(filename = "BGPC前20名基因折线图.pdf", plot = P2, device = 'pdf', width = 12, height = 36, units = 'cm')
pdf('BGPC_heatmap1.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 0,
           markGenes = markGenes,
           cluster.order = c(1:4),
           ht.col.list =list(col_range = c(-2, 0, 2),col_color = c("#16ABD0", "#FFE3DC", "#ED5A40")),
           ctAnno.col = c('#B5D465',  '#8EA0C7','#EB8F6B','#D790C1'))
dev.off()
#GOBP analysis(TOP20 based on avg_log2FC values )
library(rvcheck)
library(clusterProfiler)
library(org.Hs.eg.db)
markers.all <- Seurat::FindAllMarkers(BGPC,
                                      only.pos = TRUE,
                                      min.pct = 0.15,
                                      logfc.threshold = 0.25)
GPC <- read_excel("GPC.xlsx")
Data <- GPC$GENE 
gene.df <- bitr(Data,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)
ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoBP@result,'GPC+GO.BP.result.csv')

NRG1 <- read_excel("NRG1.xlsx")
Data <- NRG1$GENE 
gene.df <- bitr(Data,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)
ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoBP@result,'NRG1+GO.BP.result.csv')

MKI67 <- read_excel("MKI67.xlsx")
Data <- MKI67$GENE 
gene.df <- bitr(Data,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)
ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoBP@result,'MKI67+GO.BP.result.csv')

OL <- read_excel("OL.xlsx")
Data <- OL$GENE 
gene.df <- bitr(Data,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)
ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoBP@result,'OL+GO.BP.result.csv')





#mfuzz analysis
library(Mfuzz)
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
scRNA<- BGPC
scRNA$SampleID <- "SampleID"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("TOF") ] <- "TOF"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
scRNA$SampleID[ as.character(scRNA$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_4") ] <- "Con.20W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_8")] <- "Con.23W"
scRNA$SampleID[ as.character(scRNA$orig.ident) %in% c("Con_9") ] <- "Con.25W"
scRNA1 <- scRNA[,scRNA@meta.data[["SampleID"]] %in% c("TOF","Con.17W","Con.20W", "Con.21W","Con.22W","Con.23W","Con.25W")]
scRNA1=subset(scRNA1,features= rownames(scRNA1@assays$RNA@scale.data))
scRNA1@meta.data[["SampleID"]]<-factor(scRNA1@meta.data[["SampleID"]], levels=c("TOF","Con.17W","Con.20W", "Con.21W","Con.22W","Con.23W","Con.25W"))
scRNA1=ScaleData(scRNA1)
age.averages <- AverageExpression(scRNA1,group.by = "SampleID")
df1 <- age.averages[["integrated"]]
mat <- as.matrix(df1)
head(mat)
dt <- new("ExpressionSet",exprs = mat)
dim(dt)
dt.r <- filter.NA(dt, thres=0.25)
dt.f <- fill.NA(dt.r,mode="mean")
tmp <- filter.std(dt.f,min.std=0)
dt.s <- standardise(tmp)
df.s <- dt.s@assayData$exprs
m1 <- mestimate(dt.s)
set.seed(007)
cl <- mfuzz(dt.s,c=6,m=m1)
mfuzz.plot(dt.s,cl,mfrow=c(2,3), new.window= FALSE, time.labels=colnames(dt.s))
library(RColorBrewer)
mycol <- c("#E5A84B","#EE7959","#8EC9CF")
mycolor <- colorRampPalette(mycol)(20)
pdf(" integrated_TOF_expression changes trends.pdf", width = 12,height = 10)
mfuzz.plot(dt.s,cl,mfrow=c(3,3), new.window= FALSE, time.labels=colnames(dt.s), colo = mycolor)
dev.off()

#cluster5 score
retina_gene<-read_xlsx("integrted_TOF_cluster5.xlsx")
gene<-as.list(retina_gene)
AB<-AddModuleScore(scRNA1, features = gene, ctrl = 100, name = "retina")

colnames(AB@meta.data)[8]<-"retina_Score"
#####小提琴图
AB@meta.data[["celltype"]]<-factor(AB@meta.data[["celltype"]], levels=c("GPC",	"high_NRG1", "OL",	"high_MKI67"))

library(ggpubr)
color2=c(BGPC_NRG1='#8EA0C7',BGPC_APC='#B5D465',BGPC_OPC='#D790C1',BGPC_PRO='#EB8F6B')
p.AddModuleScore <- ggviolin(AB@meta.data, x = "celltype", y = "retina_Score",
                             color = "celltype",add = 'mean_sd',fill = 'celltype',
                             add.params = list(color = "black")) + 
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  scale_color_manual(values = color2) + 
  scale_fill_manual(values = color2) +
  #theme(axis.text.x.bottom = element_text(angle = 0,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')
p.AddModuleScore
ggsave('integrated_TOF_vlo_cluster5.pdf',width = 5,height = 3)

###cytotrace analysis in TOF and Con.21W
library(reticulate)
library(CytoTRACE)
reticulate::py_install(packages = c (  "scanoramaCT" )) 
BGPC.21W <- scRNA[,scRNA@meta.data[["SampleID"]] %in% c("TOF","Con.21W")]
rm(scRNA)

exp1 <- as.matrix(BGPC.21W@assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]
results <- CytoTRACE(exp1,ncores = 1)
phenot <- BGPC.21W$celltype
phenot <- as.character(phenot)
names(phenot) <- rownames(BGPC.21W@meta.data)
emb <- BGPC.21W @reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './')
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#pseudotime analysis in TOF and Con.21W
library(monocle)
library(data.table)
library(igraph)
library(dplyr)
library(heatmaply)
library(patchwork)
library(gplots)
library(ggplot2)
library(ggsci)
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
write.table(deg,file="BGPC_21Wmonocle.DEG.xls", col.names=T, row.names=F, sep="\t", quote=F)
ordergene <- row.names(deg)[order(deg$qval)][1:503]
ordergene
write.csv(ordergene, file = "BGPC_21W_ordergene.csv", append = FALSE, quote = TRUE, sep = "", row.names = TRUE, col.names = TRUE)
cds.1 <- setOrderingFilter(cds, ordergene)
table(cds.1@featureData@data[["use_for_ordering"]])
pdf("BGPC_21W_ordergenes.pdf")
plot_ordering_genes(cds.1)
dev.off()
cds.1<- reduceDimension(cds.1, max_components = 2, method = 'DDRTree')
cds.1 <- orderCells(cds.1)
cds.1 <- orderCells(cds.1,root_state = 1)
plot_cell_trajectory(cds.1, color_by = 'Pseudotime',size = 1, show_backbone = TRUE)
# change root_state
cds.1 <- orderCells(cds.1,root_state = 2)
pdf('BGPC.21W.monocle+pseudotime.pdf', width = 7, height = 7)
plot_cell_trajectory(cds.1, color_by = 'Pseudotime',size = 1, show_backbone = TRUE)
dev.off()
pdf("BGPC.21W.monocle+SampleID.pdf",width = 7,height = 7)
plot_cell_trajectory(cds.1, color_by = "SampleID",size=1,show_backbone=TRUE)
dev.off()
pdf("BGPC.21W.monocle+State.pdf",width = 7,height = 7)
plot_cell_trajectory(cds.1, color_by = "State",size=1,show_backbone=TRUE)
dev.off()

pdf("BGPC.21W.monocle+SampleID_facet.pdf",width = 7,height = 14)
plot_cell_trajectory(cds.1, color_by = "SampleID",size=1,show_backbone=TRUE)+facet_wrap("~SampleID", nrow = 2)
dev.off()
color = c("#C389BC", "#C98F8D", "#3A5A7D","#82B289")
pdf("BGPC.21W.monocle+celltype_facet.pdf",width = 7,height = 28)
plot_cell_trajectory(cds.1, color_by = "celltype",size=1,show_backbone=TRUE)+facet_wrap("~celltype", nrow = 4)+scale_colour_manual(values = color)
dev.off()


pdf('BGPC.21W.monocle+celltype.pdf', width = 7, height = 7)
plot_cell_trajectory(cds.1, color_by = 'celltype',size = 1, show_backbone = TRUE)+scale_colour_manual(values = color)
dev.off()

color = c("#C389BC", "#C98F8D", "#3A5A7D","#82B289")
p1 <- plot_complex_cell_trajectory(cds.1, x = 1, y = 2,   color_by = "celltype")+scale_colour_manual(values = color)
dev.off
ggsave("BGPC21W+celltype1.pdf", plot = p1, width = 7, height = 7)
ordergene <- BGPC_21W_ordergene
ordergene <- ordergene$x
Time_diff <-differentialGeneTest(cds.1[ordergene,], cores = 1,
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff,"BGPC.TOF.21W.Time_diff_all.csv", row.names = F)
topNGenes <- top_n(Time_diff, n = 100, desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()
pht <- plot_pseudotime_heatmap(
  cds.1[topNGenes,],
  num_clusters = 3,
  show_rownames = T,
  return_heatmap = T
)
pht
ggsave(pht, file = "time_diff_top100.pdf", width = 5, height = 20)

pdf("Time_heatmap_select.pdf", width = 5, height = 3)
marker_genes <- row.names(subset(fData(cds.1),
                                 gene_short_name %in% c("FABP7","VIM","ID2","ID4","BCAS1","SOX10","PDGFRA","ERBB3","NRG1")))
diff_test_res <- differentialGeneTest(cds.1[marker_genes,],fullModelFormulaStr ="~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
plot_pseudotime_heatmap(cds.1[sig_gene_names,], num_clusters = 1, cores = 1, show_rownames = T,hmcols = colorRampPalette(c("#83AED0","#E94A19"))(500))
dev.off()

plot_cell_trajectory(cds.1, color_by = "State")

BEAM_res <- BEAM(cds.1[ordergene,], branch_point = 1, cores = 1,progenitor_method = "duplicate") 

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res, "BEAM_res.csv", row.names = F)
library(ClusterGVis)
BEAM_gene_selected <- row.names(subset(fData(cds.1),
                                       gene_short_name %in% c("FABP7","SPARCL1","CENPE","PEA15","ADAT1",
                                                              "BCAS1","ERBB3","PDGFRA","SIRT2","MYT1L","COL9A1",
                                                              "NRG1","DLGAP2","HPSE2","SPATA17","KCNJ6","PAPPA2","PRELID2","ABL1")))
markGenes <- c("FABP7","SPARCL1","CENPE","PEA15","ADAT1",
               "BCAS1","ERBB3","PDGFRA","SIRT2","MYT1L","COL9A1",
               "NRG1","DLGAP2","HPSE2","SPATA17","KCNJ6","PAPPA2","PRELID2","ABL1")

df <- plot_genes_branched_heatmap2(cds.1[BEAM_gene_selected,],
                                   branch_point = 1,
                                   num_clusters = 3,
                                   cores = 1,
                                   use_gene_short_name = T,
                                   show_rownames = T)
str(df)
visCluster(object = df,plot.type = "heatmap")


visCluster(object = df,plot.type = "heatmap",
           pseudotime_col = c("#82B289","grey","#3A5A7D"))

pdf(file = "branch heatmap.pdf",height = 6,width = 8)
visCluster(object = df,plot.type = "both",  line.side = "right",  markGenes = markGenes,
           markGenes.side = "left", pseudotime_col = c("#82B289","grey","#3A5A7D"),
           ht.col.list =list(col_range = c(-2, 0, 2),col_color = c("#16ABD0", "#FFE3DC", "#ED5A40")),
           ctAnno.col = c("#C389BC", "#C98F8D", "#78ABDB","#82B289"))
dev.off()



pdf("BGPC.21W.keygene_celltype.pdf",width = 4,height = 4)
genes <- row.names(subset(fData(cds.1),
                          gene_short_name %in% c('ERBB3',"PDGFRA")))
plot_genes_branched_pseudotime(cds.1[genes,],
                               branch_point = 1,
                               color_by = "celltype",
                               ncol = 1)
dev.off()


library(ggpubr)
df <- pData(cds.1)
color = c("#C389BC", "#C98F8D", "#3A5A7D","#82B289")
pdf("BGPC.21W.density_celltype.pdf", width = 6, height = 4)
Clustername_color_panel <- c("GPC" = "#C389BC", "high_NRG1" = "#C98F8D", 
                             "OL" = "#3A5A7D","high_MKI67" = "#82B289")
ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) + 
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2() + 
  scale_fill_manual(name = "", values = Clustername_color_panel) + 
  scale_color_manual(name = "", values = Clustername_color_panel)
dev.off()
