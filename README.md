# HUADA-spatial-trascriptome
This script is used for recording of WZB spatial transcriptome analysis

## Sample
The samples are from hypophysoma, and one HUADA spatial transcriptomics chip contains 2 samples from patient 151 and 160 respectively.

The samples have their corresponding scRNAseq data

## Process
Start from SS200000893TL_B3.TissueCut.gem

Working Directory:
  1-2: ~/prj/shell/HUADA/WZB (PIHPC)

### 1. transvert gem to gef 
```shell
singularity exec ../SAW_v4.1.0.sif gefTools bgef -i ./SS200000893TL_B3.TissueCut.gem -o ./SS200000893TL_B3.TissueCut200.gef -b 1,200
##-b set bins for spot
```

### 2. spatial cluster by STOmics
This step will cluster spatial data by leiden, and output a .h5ad file.
``` shell
singularity exec ../SAW_v4.1.0.sif spatialCluster -i ./SS200000893TL_B3.TissueCut200.gef -o ./SS200000893TL_B3.TissueCut50.h5ad -s 50
```

### batch check 

![image](https://user-images.githubusercontent.com/49186667/218060917-0ef755d6-2ee0-4ced-b43f-b973491bf608.png)



### 3. QC and annotate clusters 

The leiden umap from step2 could be very strange (may be due to no filter step in STOmics?), present below:
<img width="591" alt="image" src="https://user-images.githubusercontent.com/49186667/217484896-6aaa1a88-bd86-4447-90eb-5f7beb7212fa.png">

So we use seurat to do filter and cluster step.

the h5ad file should be firstly transverted to h5seurat file, using SeuratDisk R package
```R
library(SeuratDisk)
Convert("SS200000893TL_B4.Tissuecut50.h5ad", dest = "h5seurat", overwrite = F)
```

#### cluster marker test

Traditional marker and SingleR will be used for cluster annotation. 

In scripts [cluster_marker_test.R], the data will be clustered into different resolustion, ranging for o.1 to 0.9 by 0.1, and markered by traditional marker signature and it own cluster differential genes. 

The filter value and PC should selected in advance

filter and PC definition. 
```
sc <- LoadH5Seurat("./SS200000893TL_B3.TissueCut10.h5seurat",meta.data = T,misc=F)
scnew <- CreateSeuratObject(counts = sc@assays$RNA,min.cells=3,min.features=0)
scnew[["percent.mt"]] <- PercentageFeatureSet(scnew, pattern = "^MT-")
#### test filter value
VlnPlot(scnew, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
scnew <- subset(scnew, subset = nFeature_RNA > 500 & nFeature_RNA < 1500 & percent.mt < 30)
scnew <- NormalizeData(scnew,assay = "RNA",normalization.method = "LogNormalize")
scnew <- FindVariableFeatures(object =scnew, selection.method = "vst",assay = "RNA",nfeatures = 3000)
scnew <- ScaleData(object = scnew,assay = "RNA")#,vars.to.regress = c('Mito.percent','Ribo.percent'))
scnew <- RunPCA(object = scnew,assay = "RNA",features = VariableFeatures(object = scnew))
#### test PC
scnew <- JackStraw(scnew, num.replicate = 100)
scnew <- ScoreJackStraw(scnew, dims = 1:20)
JackStrawPlot(scnew, dims = 1:20)
```
In this dataset, DEMO sample in bin50, should set scnew<-RunTSNE(scnew,check_duplicates = FALSE), in [cluster_marker_test.R]

The output files in this projects stored in outputs


