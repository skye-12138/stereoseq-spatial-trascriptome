# HUADA-spatial-trascriptome
This script is used for recording of WZB spatial transcriptome analysis

## Sample
The samples are from hypophysoma, and one HUADA spatial transcriptomics chip contains 2 samples from patient 151 and 160 respectively.

The samples have its corresponding scRNAseq data

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

### 3. QC and annotate clusters 
The leiden umap from step2 could be very strange (may be due to no filter step in STOmics?), present below:
<img width="591" alt="image" src="https://user-images.githubusercontent.com/49186667/217484896-6aaa1a88-bd86-4447-90eb-5f7beb7212fa.png">

So we use seurat to do filter and cluster step.

#### cluster marker test


