# Turn TF score table into a Seurat Object

countTable_bis=read.table("/Data/score_table_top50_TF.tsv", header=T, sep = "\t", row.names = 1, stringsAsFactors = F)
save(countTable_bis, file = '/Data/scoreMatrix_TF_bis.Rda')



load("/Data/scoreMatrix_TF_bis.Rda")

final_meta <- read.table("/Data/Annotated_dataset_metadata.tsv",sep = "\t",header = T)
final_meta <- final_meta %>% mutate(index=gsub("-",".",rownames(final_meta))) %>% select(index,CellType.Corrected,Position,Donor,Sample,Method,UMI.Count,Expressed.Genes,Percent.Mitochond.)
rownames(final_meta) <- final_meta$index

all(colnames(matrice)==final_meta$index)

matrice <- as(as.matrix(matrice) ,"dgCMatrix")

TF_matrix_Seurat <- CreateSeuratObject(counts = matrice, project = "TF_matrix_HCA", min.cells = 0, meta.data = final_meta)

Idents(TF_matrix_Seurat) <- "CellType.Corrected"
TF_matrix_Seurat.markers <- FindAllMarkers(object = TF_matrix_Seurat, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.5)
A <- FindMarkers(TF_matrix_Seurat,ident.1 = "Secretory N",ident.2 = "Secretory",min.pct =0.4, logfc.threshold = 0.5)
B <- FindMarkers(TF_matrix_Seurat,ident.1 = "Multiciliated N",ident.2 = "Multiciliated",min.pct =0.4, logfc.threshold = 0.5)
C <- FindMarkers(TF_matrix_Seurat,ident.1 = "Suprabasal N",ident.2 = "Suprabasal",min.pct =0.4, logfc.threshold = 0.5)


VlnPlot(TF_matrix_Seurat, features = c("SIX3"), ncol=1)  


library(ggplot2)

plots <- VlnPlot(object = TF_matrix_Seurat,cols = cell_type_color,
                 features = c('FOXA3'),
                 pt.size = 1, 
                 combine = FALSE)

for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + geom_boxplot() + theme(legend.position = 'none')
}
CombinePlots(plots)


# save(A,file="Secretory_TFTopTable.Rda")
# save(B,file="MCC_TFTopTable.Rda")
# save(C,file="Suprabasal_TFTopTable.Rda")
# save(TF_matrix_Seurat.markers,file="TFscore_TopTable.Rda")


#Embedding/Clustering sur les scores !

TF_matrix_Seurat <- ScaleData(object = TF_matrix_Seurat, features = rownames(x = TF_matrix_Seurat),verbose=F)
TF_matrix_Seurat <- RunPCA(object = TF_matrix_Seurat,features = rownames(x = TF_matrix_Seurat), verbose = FALSE)
TF_matrix_Seurat <- ProjectDim(object = TF_matrix_Seurat,verbose = 0)

DimHeatmap(object = TF_matrix_Seurat, dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(object = TF_matrix_Seurat, dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(object = TF_matrix_Seurat, dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(object = TF_matrix_Seurat, dims = 19:24, cells = 500, balanced = TRUE)

ElbowPlot(object = TF_matrix_Seurat,ndims = 30)

TF_matrix_Seurat <- RunTSNE(object = TF_matrix_Seurat, dims = 1:20)
TF_matrix_Seurat <- FindNeighbors(object = TF_matrix_Seurat, dims = 1:20,verbose = 0,k.param=10)
TF_matrix_Seurat <- FindClusters(object = TF_matrix_Seurat, resolution = 0.1,verbose = 0)
DimPlot(object = TF_matrix_Seurat, reduction = 'tsne')
TF_matrix_Seurat <- RunUMAP(object = TF_matrix_Seurat, reduction = "pca", dims = 1:20,  n.neighbors = 20)
DimPlot(object = TF_matrix_Seurat, reduction = 'umap')
p2 <-DimPlot(object = TF_matrix_Seurat, reduction = 'umap',group.by = "Donor")
p3 <-DimPlot(object = TF_matrix_Seurat, reduction = 'umap',group.by="CellType.Corrected",cols = cell_type_color)+NoLegend()
plot_grid(p2,p3)

save(TF_matrix_Seurat,file="TFscore_matrix.Rda")
