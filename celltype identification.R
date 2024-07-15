 options(stringsAsFactors = F) 
library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(Matrix)
library(gtable)
library(gridExtra)
library(patchwork)
library(EnhancedVolcano)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(fgsea)
library(enrichplot)
library(pheatmap)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(DOSE)
library(ggalluvial)
library(DoubletFinder)
library(gtable)
library(harmony)
library(clustree)
####library integration
data <- Read10X(data.dir = '43/6879-1-220919',gene.column = 1) 
FC43_1 = CreateSeuratObject(counts = data,project  = 'FC43_1', min.cells = 3, min.features = 200)  
data <- Read10X(data.dir = '43/6879-2-220919',gene.column = 1) 
FC43_2= CreateSeuratObject(counts = data,project  = 'FC43_2', min.cells = 3, min.features = 200)  
TOF = merge(FC43_1,FC43_2,add.cell.ids = c("FC43_1", "FC43_2"),project="TOF",merge.data=FALSE) # FALSE参数
intergene=intersect(rownames(FC43_1),rownames(FC43_2))
TOF=subset(TOF, features=intergene)
levels(TOF)
ids=c('TOF','TOF')
names(ids) <- levels(TOF)
TOF <- RenameIdents(TOF, ids)   
TOF$orig.ident <- Idents(TOF)
save(TOF,file="TOF.RData")

data <- Read10X(data.dir = "48/6801-1-220927",gene.column = 1)
FC48_1 = CreateSeuratObject(counts = data,project  = 'FC48_1', min.cells = 3, min.features = 200)  
data <- Read10X(data.dir = '48/6801-2-220927',gene.column = 1) 
FC48_2= CreateSeuratObject(counts = data,project  = 'FC48_2', min.cells = 3, min.features = 200)  
Con_1 = merge(FC48_1,FC48_2,add.cell.ids = c("FC48_1", "FC48_2"),project="Con_1",merge.data=FALSE) # FALSE参数
intergene=intersect(rownames(FC48_1),rownames(FC48_2))
Con_1=subset(Con_1, features=intergene)
levels(Con_1)
ids=c('Con_1','Con_1')
names(ids) <- levels(Con_1)
Con_1 <- RenameIdents(Con_1, ids)     
Con_1$orig.ident <- Idents(Con_1)
save(Con_1,file="Con_1.RData")



data <- Read10X(data.dir = "49/6791-1-220927",gene.column = 1)
FC49_1 = CreateSeuratObject(counts = data,project  = 'FC49_1', min.cells = 3, min.features = 200)  
data <- Read10X(data.dir = '49/6791-2-220927',gene.column = 1) 
FC49_2= CreateSeuratObject(counts = data,project  = 'FC49_2', min.cells = 3, min.features = 200)  
Con_2 = merge(FC49_1,FC49_2,add.cell.ids = c("FC49_1", "FC49_2"),project="Con_2",merge.data=FALSE) # FALSE参数
inter12=intersect(rownames(FC49_1),rownames(FC49_2))
Con_2 = subset(Con_2, features=inter12)
levels(Con_2)
ids=c('Con_2','Con_2')
names(ids) <- levels(Con_2)
Con_2 <- RenameIdents(Con_2, ids)      
Con_2$orig.ident <- Idents(Con_2)
save(Con_2,file="Con_2.RData")

data <- Read10X(data.dir = "50/6799-1-220927",gene.column = 1)
FC50_1 = CreateSeuratObject(counts = data,project  = 'FC50_1', min.cells = 3, min.features = 200)  
data <- Read10X(data.dir = '50/6799-2-220927',gene.column = 1) 
FC50_2= CreateSeuratObject(counts = data,project  = 'FC50_2', min.cells = 3, min.features = 200)  
Con_3 = merge(FC50_1,FC50_2,add.cell.ids = c("FC50_1", "FC50_2"),project="Con_3",merge.data=FALSE) # FALSE参数
inter12=intersect(rownames(FC50_1),rownames(FC50_2))
Con_3 = subset(Con_3, features=inter12)
levels(Con_3)
ids=c('Con_3','Con_3')
names(ids) <- levels(Con_3)
Con_3 <- RenameIdents(Con_3, ids)     
Con_3$orig.ident <- Idents(Con_3)
save(Con_3,file="Con_3.RData")

data <- Read10X(data.dir = "53/6794-1-220928",gene.column = 1)
FC53_1 = CreateSeuratObject(counts = data,project  = 'FC53_1', min.cells = 3, min.features = 200)  
data <- Read10X(data.dir = '53/6794-2-220928',gene.column = 1) 
FC53_2= CreateSeuratObject(counts = data,project  = 'FC53_2', min.cells = 3, min.features = 200)  
Con_4 = merge(FC53_1,FC53_2,add.cell.ids = c("FC53_1", "FC53_2"),project="Con_4",merge.data=FALSE) # FALSE参数
inter12=intersect(rownames(FC53_1),rownames(FC53_2))
Con_4 = subset(Con_4, features=inter12)
levels(Con_4)
ids=c('Con_4','Con_4')
names(ids) <- levels(Con_4)
Con_4 <- RenameIdents(Con_4, ids)     
Con_4$orig.ident <- Idents(Con_4)
save(Con_4,file="Con_4.RData")

data <- Read10X(data.dir = "51/6783-1-220919",gene.column = 1)
FC51_1 = CreateSeuratObject(counts = data,project  = 'FC51_1', min.cells = 3, min.features = 200)  
data <- Read10X(data.dir = '51/6783-2-220919',gene.column = 1) 
FC51_2= CreateSeuratObject(counts = data,project  = 'FC51_2', min.cells = 3, min.features = 200)  
Con_5 = merge(FC51_1,FC51_2,add.cell.ids = c("FC51_1", "FC51_2"),project="Con_5",merge.data=FALSE) # FALSE参数
inter12=intersect(rownames(FC51_1),rownames(FC51_2))
Con_5 = subset(Con_5, features=inter12)
levels(Con_5)
ids=c('Con_5','Con_5')
names(ids) <- levels(Con_5)
Con_5 <- RenameIdents(Con_5, ids)    
Con_5$orig.ident <- Idents(Con_5)
save(Con_5,file="Con_5.RData")


data <- Read10X(data.dir = "77/6894-1-220915",gene.column = 1)
FC77_1 = CreateSeuratObject(counts = data,project  = 'FC77_1', min.cells = 3, min.features = 200)  
data <- Read10X(data.dir = '77/6894-2-220915',gene.column = 1) 
FC77_2= CreateSeuratObject(counts = data,project  = 'FC77_2', min.cells = 3, min.features = 200)  
Con_6 = merge(FC77_1,FC77_2,add.cell.ids = c("FC77_1", "FC77_2"),project="Con_6",merge.data=FALSE) # FALSE参数
inter12=intersect(rownames(FC77_1),rownames(FC77_2))
Con_6 = subset(Con_6, features=inter12)
levels(Con_6)
ids=c('Con_6','Con_6')
names(ids) <- levels(Con_6)
Con_6 <- RenameIdents(Con_6, ids)      
Con_6$orig.ident <- Idents(Con_6)
save(Con_6,file="Con_6.RData")

data <- Read10X(data.dir = "78/6895-1-220915",gene.column = 1)
FC78_1 = CreateSeuratObject(counts = data,project  = 'FC78_1', min.cells = 3, min.features = 200)  
data <- Read10X(data.dir = '78/6895-2-220915',gene.column = 1) 
FC78_2= CreateSeuratObject(counts = data,project  = 'FC78_2', min.cells = 3, min.features = 200)  
Con_7 = merge(FC78_1,FC78_2,add.cell.ids = c("FC78_1", "FC78_2"),project="Con_7",merge.data=FALSE) # FALSE参数
inter12=intersect(rownames(FC78_1),rownames(FC78_2))
Con_7= subset(Con_7, features=inter12)
levels(Con_7)
ids=c('Con_7','Con_7')
names(ids) <- levels(Con_7)
Con_7 <- RenameIdents(Con_7, ids)       
Con_7$orig.ident <- Idents(Con_7)
save(Con_7,file="Con_7.RData")

data <- Read10X(data.dir = "X7/4412-1-230718",gene.column = 1)
X7_1 = CreateSeuratObject(counts = data,project  = 'X7_1', min.cells = 3, min.features = 200)  
data <- Read10X(data.dir = 'X7/4412-2-230718',gene.column = 1) 
X7_2= CreateSeuratObject(counts = data,project  = 'X7_2', min.cells = 3, min.features = 200)  
Con_8 = merge(X7_1,X7_2,add.cell.ids = c("X7_1", "X7_2"),project="Con_8",merge.data=FALSE) # FALSE参数
inter12=intersect(rownames(X7_1),rownames(X7_2))
Con_8= subset(Con_8, features=inter12)
levels(Con_8)
ids=c('Con_8','Con_8')
names(ids) <- levels(Con_8)
Con_8 <- RenameIdents(Con_8, ids)       
Con_8$orig.ident <- Idents(Con_8)
save(Con_8,file="Con_8.RData")

data <- Read10X(data.dir = "X8/4413-1-230718",gene.column = 1)
X8_1 = CreateSeuratObject(counts = data,project  = 'X8_1', min.cells = 3, min.features = 200)  
data <- Read10X(data.dir = 'X8/4413-2-230718',gene.column = 1) 
X8_2= CreateSeuratObject(counts = data,project  = 'X8_2', min.cells = 3, min.features = 200)  
Con_9 = merge(X8_1,X8_2,add.cell.ids = c("X8_1", "X8_2"),project="Con_9",merge.data=FALSE) # FALSE参数
inter12=intersect(rownames(FC78_1),rownames(FC78_2))
Con_9= subset(Con_9, features=inter12)
levels(Con_9)
ids=c('Con_9','Con_9')
names(ids) <- levels(Con_9)
Con_9 <- RenameIdents(Con_9, ids)       
Con_9$orig.ident <- Idents(Con_9)
save(Con_9,file="Con_9.RData")

obj.list = list()

obj.list[['TOF']]= TOF
obj.list[['Con_1']]= Con_1
obj.list[['Con_2']]= Con_2
obj.list[['Con_3']]= Con_3
obj.list[['Con_4']]= Con_4
obj.list[['Con_5']]= Con_5
obj.list[['Con_6']]= Con_6
obj.list[['Con_7']]= Con_7
obj.list[['Con_8']]= Con_8
obj.list[['Con_9']]= Con_9


obj.list <- lapply(obj.list,function(x) {
  NormalizeData(x)
})
obj.list <- lapply(obj.list, function(x) {
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
obj.list
obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
obj.combined <- IntegrateData(anchorset = obj.anchors, dims = 1:30)
DefaultAssay(obj.combined) <- "integrated"
obj.combined <- ScaleData(obj.combined, features = rownames(obj.combined))
obj.combined <- RunPCA(obj.combined, features = VariableFeatures(object = obj.combined))
print(obj.combined[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(obj.combined, dims = 1:5, reduction = "pca",ncol=2)
DimPlot(obj.combined, reduction = "pca")
DimHeatmap(obj.combined, dims = 1:15, cells = 500, balanced = TRUE)

obj.combined <- JackStraw(obj.combined, num.replicate = 100)
obj.combined <- ScoreJackStraw(obj.combined, dims = 1:20)
plot1 <- JackStrawPlot(obj.combined, dims = 1:20)
plot2 <- ElbowPlot(obj.combined,ndims = 50)

dimNums = 20 
obj.combined <- RunUMAP(obj.combined, dims = 1:dimNums)
obj.combined <- FindClusters(obj.combined,resolution = 0.6)
Idents(obj.combined) <-  "integrated_snn_res.0.6"
obj.combined$seurat_clusters <- obj.combined@active.ident
doc_makers=list(
  'ExN'=c('RBFOX1','SATB2','SLC17A7',"CAMK2A",'NEUROD2','NEUROD6','CNTN5','RELN'),
  'InN'=c('GAD1','GAD2', 'SLC3AA1'),
  'Astro'=c('AQP4','GFAP','SLC1A3','SLC1A2','SOX9'),
  'Micro'=c('INPP5D','CSF1R','LGMN','P2RY12','IFNGR1','C1QA','CXCR1','PTRRC'),
  'Oligo'=c('PLP1','MOG','MOBP','MBP','BCAS1','NFASC','TCF7L2','SIRT2'),
  'OPC'=c('PDGFRA','PCDH15','NEU4','MYT1','VCAN','SOX10', 'VCAN','OLG1', 'OLG2'),
  'Endo'=c('VWF','CLDN5','FLT1'),
  'NPCS'=c('PAX6','SFRP1','DCX')
)
genes_to_check=unique(unlist(doc_makers))
DotPlot(obj.combined, assay = "RNA", features = genes_to_check, cols = c("blue", "red"), col.min = 0,col.max = 2) + coord_flip()+ RotatedAxis()



obj.combined2 <- obj.combined[,obj.combined@meta.data[["seurat_clusters"]] %in% c("1",	"2",	"3", "4",	"5",	"10",	"11")]
obj.combined2<- RunUMAP(obj.combined2, dims = 1:20)
p8 <- DimPlot(obj.combined2, reduction = "umap", group.by = "ident", pt.size=0.1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
doc_makers=list(
  'ExN'=c('RBFOX1','SATB2','SLC17A7',"CAMK2A",'NEUROD2','NEUROD6','CNTN5','RELN'),
  'InN'=c('GAD1','GAD2', 'SLC3AA1'),
  'Astro'=c('AQP4','GFAP','SLC1A3','SLC1A2','SOX9'),
  'Micro'=c('INPP5D','CSF1R','LGMN','P2RY12','IFNGR1','C1QA','CXCR1','PTRRC'),
  'Oligo'=c('PLP1','MOG','MOBP','MBP','BCAS1','NFASC','TCF7L2','SIRT2'),
  'OPC'=c('PDGFRA','PCDH15','NEU4','MYT1','VCAN','SOX10', 'VCAN','OLG1', 'OLG2'),
  'Endo'=c('VWF','CLDN5','FLT1'),
  'NPCS'=c('PAX6','SFRP1','DCX')
)
genes_to_check=unique(unlist(doc_makers))
DotPlot(obj.combined2, assay = "RNA", features = genes_to_check, cols = c("blue", "red"), col.min = -1,col.max = 2) + coord_flip()+ RotatedAxis()
###celltype identification
current.cluster.ids <- c("1",	"10",	"11", "2",	"3",	"4",	"5")
new.cluster.ids <- c("ExN","GC",	"GC",	"ExN",	"ExN",	"InN",	"ExN") 
obj.combined2@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(obj.combined2@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(obj.combined2@meta.data$celltype)
current.cluster.ids <- c("1",	"10",	"11", "2",	"3",	"4",	"5")
names(new.cluster.ids) <- levels(obj.combined2)
obj.combined2 <- RenameIdents(obj.combined2, new.cluster.ids)
DimPlot(obj.combined2, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
save(obj.combined2,file="obj.combined2.RData")

p9 <- DimPlot(obj.combined2, cols = c('#B5D465',  '#8EA0C7', '#D790C1'), reduction = "umap", group.by = "ident", pt.size=0.1, label = T,label.size=5, repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "UMAP+celltype.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')

obj.combined2$group=str_replace(obj.combined2$orig.ident,"_//d","")
Idents(obj.combined2)='group'
obj.combined2$group <- factor(as.character(obj.combined2$group), levels=c("TOF", "Con_1", "Con_2", "Con_3", "Con_4", "Con_5" ,"Con_6" ,"Con_7", "Con_8" ,"Con_9"))
Idents(obj.combined2) <- obj.combined2$group
p1 <- VlnPlot(obj.combined2, group.by = "group",features = c("nFeature_RNA"), cols = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",  "#D9E6B7","#BC80BD"), pt.size =0,ncol=1)
ggsave("ALL_control_nFeature_RNA.pdf",  plot = p1,width = 20, height = 8, units = "cm")
p1 <- VlnPlot(obj.combined2, group.by = "group",features = c("nCount_RNA"), cols = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",  "#D9E6B7","#BC80BD"), pt.size = 0 ,ncol =1)
ggsave("ALL_control_nCount_RNA.pdf", plot = p1,width = 20, height = 8, units = "cm")


obj.combined2$SampleID <- "SampleID"
obj.combined2$SampleID[ as.character(obj.combined2$orig.ident) %in% c("TOF") ] <- "TOF"
obj.combined2$SampleID[ as.character(obj.combined2$orig.ident) %in% c("Con_1","Con_2") ] <- "Con.21W"
obj.combined2$SampleID[ as.character(obj.combined2$orig.ident) %in% c("Con_3","Con_5") ] <- "Con.22W"
obj.combined2$SampleID[ as.character(obj.combined2$orig.ident)%in% c("Con_6","Con_7") ] <- "Con.17W"
obj.combined2$SampleID[ as.character(obj.combined2$orig.ident) %in% c("Con_4") ] <- "Con.20W"
obj.combined2$SampleID[ as.character(obj.combined2$orig.ident) %in% c("Con_8")] <- "Con.23W"
obj.combined2$SampleID[ as.character(obj.combined2$orig.ident) %in% c("Con_9") ] <- "Con.25W"
table(obj.combined2@meta.data$SampleID)
p9 <- DimPlot(obj.combined2, cols = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",  "#D9E6B7","#BC80BD"), reduction = "umap", group.by = "group", pt.size=0.01, label = F,label.size=5, repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "UMAP+group.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')
p9 <- DimPlot(obj.combined2, cols = c("#8DD3C7", "#FDB462", "#BEBADA", "#FB8072", "#80B1D3","#FFFFB3", "#B3DE69"), reduction = "umap", group.by = "SampleID",pt.size=0.01, label = F,label.size=5, repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap+SampleID.pdf", plot = p9, device = 'pdf', width = 20, height = 20, units = 'cm')
AB <- obj.combined2[,obj.combined2@meta.data[["SampleID"]] %in% c("TOF")]
p10 <- DimPlot(AB, cols = c('#D9E6B7',  '#8EA0C7', '#D790C1'), reduction = "umap", group.by = "celltype", pt.size=0.5, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "UMAP+TOF.pdf", plot = p10, device = 'pdf', width = 15, height = 15, units = 'cm')
AB <- obj.combined2[,obj.combined2@meta.data[["SampleID"]] %in% c("Con.21W")]
p11 <- DimPlot(AB, cols = c('#D9E6B7',  '#8EA0C7', '#D790C1'), reduction = "umap", group.by = "celltype", pt.size=0.5, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "UMAP+Con.21w.pdf", plot = p11, device = 'pdf', width = 15, height = 15, units = 'cm')

genes_to_check = c('RBFOX1','SATB2','NEUROD2','SLC1A3','SLC1A2', 'INPP5D','TCF7L2','SIRT2','PDGFRA','PCDH15','GAD1','ERBB4' )
library(stringr)  
genes_to_check=str_to_upper(unique(genes_to_check))
genes_to_check

th=theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
p_all_markers <- DotPlot(obj.combined2, features = genes_to_check,
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
ggsave(file="maker_exp_dotplot_obj.combined2.pdf",plot = p1,height= 6,width = 5)


###find DEGs 
DefaultAssay(obj.combined2) <- "RNA"
table(obj.combined2@meta.data$orig.ident)
obj.combined2$group=str_replace(obj.combined2$orig.ident,"_//d","")
table(obj.combined2@meta.data$group)
table(obj.combined2@meta.data$celltype)
Idents(obj.combined2)='group'
levels(obj.combined2) 
group2 = c('TOF',"Con_21","Con_21","Con_22","Con_20", "Con_22", 'Con_17', 'Con_17', 'Con_23','Con_25')
names(group2) = levels(obj.combined2)
obj.combined2 <- RenameIdents(obj.combined2, group2)
levels(obj.combined2)
obj.combined2$group2 <- Idents(obj.combined2)
obj.combined2$celltype.group2 <- paste(obj.combined2$celltype,obj.combined2$group2, sep = "_") 
Idents(obj.combined2) <- "celltype.group2"
table(obj.combined2@active.ident)
ExN.TOF_CON_21.diff <- FindMarkers(obj.combined2, ident.1 = "ExN_TOF", ident.2 = "ExN_Con_21", logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
library(tibble)
ExN.TOF_CON_21.diff<- tibble::rownames_to_column(ExN.TOF_CON_21.diff, "gene")
colnames(ExN.TOF_CON_21.diff)[3] <-"log2FoldChange"
colnames(ExN.TOF_CON_21.diff)[6] <-"padj"
write.table(ExN.TOF_CON_21.diff,file="ExN_TOF_CON21W.diff_genes.txt",quote=F,sep="\t",row.names=F,col.names=T)
InN.TOF_CON_21.diff <- FindMarkers(obj.combined2, ident.1 = "InN_TOF", ident.2 = "InN_TOF", logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
InN.TOF_CON_21.diff<- tibble::rownames_to_column(InN.TOF_CON_21.diff, "gene")
colnames(InN.TOF_CON_21.diff)[3] <-"log2FoldChange"
colnames(InN.TOF_CON_21.diff)[6] <-"padj"
write.table(InN.TOF_CON_21.diff,file="InN.TOF_CON_21.diff_genes.txt",quote=F,sep="\t",row.names=F,col.names=T)
GC.TOF_CON_21.diff <- FindMarkers(obj.combined2, ident.1 = "GC_TOF", ident.2 = "GC_Con_21", logfc.threshold = 0.25, min.pct = 0.15, only.pos = F, test.use = "wilcox")
GC.TOF_CON_21.diff<- tibble::rownames_to_column(GC.TOF_CON_21.diff, "gene")
colnames(GC.TOF_CON_21.diff)[3] <-"log2FoldChange"
colnames(GC.TOF_CON_21.diff)[6] <-"padj"
write.table(GC.TOF_CON_21.diff,file="GC.TOF_CON_21.diff_genes.txt",quote=F,sep="\t",row.names=F,col.names=T)

if(!require(rvcheck))devtools::install_version("rvcheck", version = "0.1.8", repos = "http://cran.us.r-project.org")
if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
if(!require(org.Hs.eg.db))BiocManager::install("org.Hs.eg.db")
if(!require(rvcheck))BiocManager::install("rvcheck")
library(rvcheck)
library(clusterProfiler)
library(org.Hs.eg.db)
if(!require(clusterProfiler))BiocManager::install("clusterProfiler")
if(!require(igraph))BiocManager::install("igraph")
installed.packages("igraph")
###GOBP analysis
Data <- read.csv('ExN_DOWN.txt',sep= '')
Data <- ExN下调
mygene <- Data$gene 
gene.df <- bitr(mygene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)

ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoBP@result,'EXN_DOWNGO.BP.result.csv')

Data <- read.csv('GC_DOWN.txt',sep= '')

mygene <- Data$gene
gene.df <- bitr(mygene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)

ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoBP@result,'GC_DOWNGO.BP.result.csv')

Data <- read.csv('InN_DOWN.txt',sep= '')

mygene <- Data$gene 
gene.df <- bitr(mygene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)

ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoBP@result,'InN_DOWNGO.BP.result.csv')

##pseudotime analysis
scRNA<- obj.combined2
scRNA[['integrated_snn_res.0.6']] <- NULL
head(scRNA@meta.data)
scRNA$celltype_SampleID <- paste0(scRNA$celltype,"*",scRNA$SampleID)
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
save(cds,file = 'input_cds.Rdata')
diff <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr='~celltype', cores=1)
head(diff)
deg <- subset(diff, qval < 0.05)
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
write.table(deg,file="monocle.DEG.xls", col.names=T, row.names=F, sep="\t", quote=F)
ordergene <- row.names(deg)[order(deg$qval)][1:2000]
ordergene
cds.1<-cds
write.csv(ordergene, file = "ordergene.csv", append = FALSE, quote = TRUE, sep = "", row.names = TRUE, col.names = TRUE)
cds <- setOrderingFilter(cds, ordergene)
table(cds@featureData@data[["use_for_ordering"]])
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = 1)
pdf('monocle+pseudotime.pdf', width = 7, height = 7)
plot_cell_trajectory(cds, color_by = 'Pseudotime',size = 1, show_backbone = TRUE)
dev.off()
color = c('#D9E6B7',  '#8EA0C7', '#D790C1')
pdf('monocle+celltype.pdf', width = 7, height = 7)
plot_cell_trajectory(cds, color_by = 'celltype',size = 1, show_backbone = TRUE)+
  scale_colour_manual(
    values = color)
dev.off()
color = c("#F2A49A","#80B9C8")
pdf('monocle+SampleID.pdf', width = 7, height = 7)
plot_cell_trajectory(cds, color_by = 'SampleID',size = 1, show_backbone = TRUE)+
  scale_colour_manual(
    values = color)
dev.off()
pdf(file = "monocle+树形图+ExN+GC+InN.pdf", width = 10, height = 10)
plot_complex_cell_trajectory(cds, x = 1, y = 2,   color_by = "SampleID")
dev.off()