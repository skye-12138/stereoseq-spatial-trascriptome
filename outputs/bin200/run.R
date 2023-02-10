library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(patchwork)
library(SingleR)
library(celldex)
library(cowplot)
library(paletteer)
######
###matrix for marker usage should be generated from seurat FindAllMarkers
##matrix<-FindAllMarkers(seurat_obj)
statistic_plot<-function(matrix,seurat_obj=NULL,type="",usage="",condition=NULL,save=FALSE,top=10){
  if(usage == "cellratio"){
    if(type == "bar"){
      P<-ggplot(matrix) +
        aes(x = get(sample), y = Percentage, fill = cellType) +
        geom_col() +
        scale_fill_hue(direction = 1) +
        theme_minimal()+
        labs(x=sample)+
        scale_fill_manual(values = paletteer::paletteer_d("ggsci::category20_d3",n=length(unique(seurat_obj@meta.data[[celltype]]))))
      if(is.null(condition) == F){
        P<-P+facet_wrap(.~condition,scales = 'free')
      }
    }
    else if(type == "pie"){
      P<-ggplot(matrix) +
        aes(x= sub, fill = cellType, weight = num) +
        geom_bar(position = "fill") +
        scale_fill_hue() +coord_polar(theta = 'y')+
        theme_pubclean()+
        theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y = element_blank(),axis.text.y = element_blank()) +
        scale_fill_manual(values = paletteer_d("ggsci::category20_d3",n=length(unique(matrix$cellType))),labels = paste(unique(matrix$cellType),"(",scales::percent(matrix$ratio),")",sep = ""))
    }
  }######plot for marker presentation####Seurat object is necessary in this section
  else if(usage == "marker"){
    matrix %>%
      group_by(cluster) %>%
      top_n(n = top, wt = avg_log2FC) -> topmarker
    if(type == "heatmap"){
      P<-DoHeatmap(seurat_obj,assay = "RNA",slot = "scale.data",features=as.character(unique(topmarker$gene)),group.by = "cellType",group.colors = paletteer_d("ggsci::category20_d3",n=length(levels(Idents(seurat_obj)))))+scale_fill_gradientn(colors=c("#26456E","white","#C53E01"))
    }else if(type == "heatmap_av"){
      avg<-AverageExpression(seurat_obj,assay = "RNA",features = as.character(unique(topmarker$gene)),group.by = "cellType",slot = "scale.data")
      avg<-avg$RNA
      P<-pheatmap::pheatmap(avg,scale = "row",cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c("#1A237E","white","#ff7f00"))(length(seq(0,2,by = 0.2))),breaks =seq(0,2,by = 0.2),border=FALSE)
    }else if(type == "vln"){
      P<-VlnPlot(seurat_obj, features = as.character(unique(topmarker$gene)), stack = TRUE,cols = rep(paletteer_d("ggsci::category20_d3",n=length(levels(Idents(seurat_obj)))),each=top)) +
        theme(legend.position = "none")
    }else if(type == "dot"){
      P<-DotPlot(seurat_obj, features = unique(topmarker$gene),
                 assay='RNA' )
    }
  }
  ###can add other usage and plot script here 
  return(P)
}
######
plot_feature<-function(seu,feature,name){
  seu<-AddModuleScore(seu,features = feature,name=names(feature))
  colnames(seu@meta.data)[c(length(colnames(seu@meta.data)) - length(names(feature)) + 1):length(colnames(seu@meta.data))]<-names(feature)
  P1<-FeaturePlot(seu,features = names(feature),min.cutoff = "q10",reduction = "spatial")
  P2<-FeaturePlot(seu,features = names(feature),min.cutoff = "q10",reduction = "tsne")
  png(paste(name,".png",sep=""),width = 2200,height = 1000)
  print(cowplot::plot_grid(plotlist = list(P1,P2)))
  dev.off()
  ####group.by="leiden" only use in STOmics situation
  P3<-DotPlot(seu,features = names(feature),cols = c("dark blue","red"))+theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))+ggtitle(label = name)
  return(P3)
}

######
########
########load data
#######
hpca.se <- HumanPrimaryCellAtlasData()
Blue.se=BlueprintEncodeData()
sc <- LoadH5Seurat("./SS200000893TL_B3.TissueCut200.h5seurat",meta.data = T,misc=F)
scnew <- CreateSeuratObject(counts = sc@assays$RNA,min.cells=3,min.features=0)
scnew@reductions$spatial<-sc@reductions$spatial
scnew[["percent.mt"]] <- PercentageFeatureSet(scnew, pattern = "^MT-")
VlnPlot(scnew, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
scnew<-subset(scnew, subset = nFeature_RNA > 4000 & nFeature_RNA < 7000 & percent.mt < 25)
scnew <- NormalizeData(scnew,assay = "RNA",normalization.method = "LogNormalize")
scnew <- FindVariableFeatures(object =scnew, selection.method = "vst",assay = "RNA",nfeatures = 3000)
scnew <- ScaleData(object = scnew,assay = "RNA")
scnew <- RunPCA(object = scnew,assay = "RNA",features = VariableFeatures(object = scnew))
scnew <- FindNeighbors(scnew, reduction = "pca", dims = 1:17)
for(i in seq(0.1,0.9,by=0.1)){
  scnew<-FindClusters(scnew,algorithm = 1,resolution = i) 
}
scnew <- RunUMAP(object = scnew, dims = 1:17)
scnew<-RunTSNE(scnew)
####marker plot
marker_sig<-list(Tumor_POU1F1=c("POU1F1","GH1","PRL","PRLR"),
                 Tumor_TBX19=c("TBX19","GZMK","MAPK12","AVPR1B","CXCL13","AKT3","POMC"),
                 Tumor_NR5A1=c("NR5A1","EGR4","EGR2","FSHB","IDH1"),
                 Marco_DC=c("C1QA","C1QB","C1QC","CCL3","MS4A7","CCL3L1","CCL4L2"),
                 Neutrophil=c("S100A8","S100A9","FCGR3B","CSF3R","CXCR2"),
                 Bcells=c("MS4A1","CD79A","IGHM","IGKC"),
                 T_NKcells=c("CD3E","NKG7","IL2RB","GZMA","CCL5"),
                 MAST=c("KIT","MS4A2","CPA3","TPSB2"),
                 Endothelial=c("VWF","PLVAP","STC1"),
                 Epithelial=c("CXCL17","AGR2","SLPI"),
                 Fibroblast=c("DCN","FN1","MYH11","COL1A2")
)
for(i in names(marker_sig)){
  Pfeature1<-FeaturePlot(scnew,features = marker_sig[[i]],reduction = "tsne",repel = T)+ggtitle(label = i)
  Pfeature2<-FeaturePlot(scnew,features = marker_sig[[i]],reduction = "spatial",repel = T)+ggtitle(label = i)
  Pfeature<-plot_grid(plotlist = list(Pfeature1,Pfeature2),nrow = 2)
  png(paste(i,"_subfeaturemarker.png",sep=""),height = 1000,width = 2200)
  print(Pfeature)
  dev.off()
}
####singleR
test_for_singler<-GetAssayData(scnew, slot="data")
singler1<-SingleR(test= test_for_singler , ref = list(hp=hpca.se,bl=Blue.se) ,labels = list(hpca.se$label.main,Blue.se$label.main))
scnew$single_inside_nocluster<-singler1$labels
####
for(i in seq(0.1,0.9,by=0.1)){
  Idents(scnew)<-scnew[[paste("RNA_snn_res.",i,sep = "")]]
  pmarker_sig<-plot_feature(scnew,marker_sig,"marker_sig")
  png(paste(i,"marker_sig.png",sep="_"),width = 1500,height = 1500)
  print(pmarker_sig)
  dev.off()
  marker<-FindAllMarkers(scnew)
  write.table(marker,paste(i,"marker.txt",sep="_"),quote = F)
  scnew$cellType<-scnew[[paste("RNA_snn_res.",i,sep = "")]]
  P<-statistic_plot(matrix = marker,seurat_obj = scnew,type="heatmap",usage="marker")
  pdf(paste(i,"heatmap.pdf",sep = "_"))
  print(P)
  dev.off()
  ####singleR-singlecell
  P1<-DimPlot(scnew, group.by=paste("RNA_snn_res.",i,sep = ""), label=T, label.size=5, reduction='tsne')
  P2<-DimPlot(scnew, group.by=paste("RNA_snn_res.",i,sep = ""), label=T, label.size=5, reduction='umap')
  P3<-DimPlot(scnew, group.by="single_inside_nocluster", label=T, label.size=5, reduction='tsne')
  P4<-DimPlot(scnew, group.by="single_inside_nocluster", label=T, label.size=5, reduction='umap')
  ####singleR-cluster
  clusters=scnew@meta.data[[paste("RNA_snn_res.",i,sep = "")]]
  singler2<-SingleR(test= test_for_singler,ref = list(hp=hpca.se,bl=Blue.se) ,labels = list(hpca.se$label.main,Blue.se$label.main),clusters=clusters)
  celltype = data.frame(ClusterID=rownames(singler2), celltype=singler2$labels, stringsAsFactors = FALSE)
  scnew@meta.data[[paste("celltype",i,sep=".")]] = "NA"
  for(n in 1:nrow(celltype)){
    scnew@meta.data[which(scnew@meta.data[[paste("RNA_snn_res.",i,sep = "")]] == celltype$ClusterID[n]),paste("celltype",i,sep=".")] <- celltype$celltype[n]}
  p1 = DimPlot(scnew, group.by=paste("celltype",i,sep="."), label=T, label.size=5, reduction='tsne')
  # umap
  p2 = DimPlot(scnew, group.by=paste("celltype",i,sep="."), label=T, label.size=5, reduction='umap')
  # 合并
  p3 = plotc <- P1+P2+P3+P4+p1+p2+patchwork::plot_layout(guides = "collect",ncol = 2)
  ggsave(paste(i,"celltype_cluster.pdf",sep="_"), p3, width=15 ,height=8)
}
saveRDS(scnew,"scnew.rds")


