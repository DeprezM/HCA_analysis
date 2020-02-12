## Fig Robustness
library(ggplot2)
library(pheatmap)
library(tidyr)

setwd("/HCA_analysis/6_Figures/Heatmaps_robust")

##########################################################################################
#                       Load datasets & Create variables                                 #
##########################################################################################
cell_type_order <- c('Cycling Basal',
                     'Basal',
                     'Suprabasal', 'Suprabasal N',
                     'Secretory', 'Secretory N',
                     'Deuterosomal', 
                     'Multiciliated', 'Multiciliated N',
                     'SMG Goblet', 'Serous', 'Rare cells', 
                     'AT1', 'AT2', 'Brush cells', 'PNEC', 'Precursor', 'Ionocyte',
                     'Endothelial', 'Fibroblast', 'Pericyte', 'Smooth muscle',
                     'B cells', 'Plasma cells', 'LT/NK', 'Mast cells', 'Dendritic', 'Monocyte', 'Macrophage')

method_color <-  c('Brushing' = '#009933',
                   'Biopsy' = '#336699')

position_color <- c('Nasal' = '#ffd966',
                    'Proximal' = '#e69138',
                    'Intermediate' = '#e06666',
                    'Distal' = '#85200c')

flexible_normalization <- function(data_in, by_row=TRUE){
  if(by_row){
    row_mean <- apply(data_in,1,mean)
    row_sd   <- apply(data_in,1,sd)
    output <- data_in
    for(i in 1:dim(data_in)[1]){
      output[i,] <- (data_in[i,] - row_mean[i])/row_sd[i]
    }
  }
  #### if by column
  if(!by_row){
    col_mean <- apply(data_in,2,mean)
    col_sd   <- apply(data_in,2,sd)
    output <- data_in
    for(i in 1:dim(data_in)[2]){
      output[,i] <- (data_in[,i] - col_mean[i])/col_sd[i]
    }
  }
  return(output)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

## Load annotated dataset ----------------------------------------------------------------
metadata <- read.table(file = "/Data/Annotated_dataset_metadata.tsv", 
                       sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
coord <- read.table(file = "/Data/Annotated_dataset_umap_coords.tsv", 
                    sep = "\t", header = T)
cluster_color <- read.csv(file = "/Data/clusterColor.tsv",
                          sep = "\t", stringsAsFactors = F, header = T)
cell_type_color <- cluster_color$color
names(cell_type_color) <- cluster_color$Cluster

metadata$x = coord$x
metadata$y = coord$y
metadata$Position <- factor(metadata$Position, levels = names(position_color))
metadata$Method <- factor(metadata$Method, levels = names(method_color))
metadata$SampleType = paste0(metadata$Position, '\n', metadata$Method)
metadata$SampleType = factor(metadata$SampleType, levels = sampletype_order)


# Residuals <100 cells from Nasal or Bronchial position labelled into another cluster
# metadata[metadata$CellType.Corrected == "Suprabasal N" & metadata$Position != "Nasal",]$CellType.Corrected <- "Suprabasal"
# metadata[metadata$CellType.Corrected == "Secretory N" & metadata$Position != "Nasal",]$CellType.Corrected <- "Secretory"
# metadata[metadata$CellType.Corrected == "Multiciliated N" & metadata$Position != "Nasal",]$CellType.Corrected <- "Multiciliated"


### Set up necessary variables --------------------------------------------------------
library("readxl")
marker_genes <- read_excel("/Data/MarkerGenesHuman.xlsx", sheet = 1)

load('/Data/SoupX_norm_dataset.Rda')

epithelial <- c("Cycling Basal", 'Basal', 'Suprabasal', 'Suprabasal N', 'Club','Secretory N', 'Deuterosomal', 
                'Multiciliated', 'Multiciliated N',
                'SMG Goblet', 'Serous', 
                'AT1', 'AT2', 'Tuft Cell / Brush cell', 'Neuro-endocrine', 'Ionocyte')
immune <- c( 'LB', 'Plasmocyte', 'Lymphocyte', 'Mastocyte', 'Dendritic Cells', 'Monocyte/Neutro', 'Macrophage')
stromal <- c('Endothelial', 'Fibroblast', 'Smooth muscle')


##########################################################################################
#                                     Heatmaps                                           #
##########################################################################################


#### EPITHELIAL -------------------------------------------------------------------------

markers <- c()
for (i in epithelial) {
  markers <- c(markers, unlist(marker_genes[,i]))
}

markers <- markers[!is.na(markers)]
markers <- markers[markers %in% rownames(countTable)]
markers_flt_count <- countTable[markers, ]

cell_type_names <- c("Cycling Basal", "Basal", "Suprabasal", "Suprabasal N", "Secretory", 
                     "Secretory N", "Deuterosomal", "Multiciliated", "Multiciliated N",
                     "SMG Goblet", "Serous",  "AT1", "AT2",  'Brush cells', 'PNEC', 'Ionocyte')
sample_names <- unique(metadata[order(metadata$Position, metadata$Method), "Sample"])
markers_count <- matrix(0, nrow = length(markers), ncol = length(cell_type_names)*length(unique(metadata$Sample)))

rownames(markers_count) <- markers
colnames_counts <- c()
i=1
for (j in 1:length(cell_type_names)){
  print(cell_type_names[j])
  for (k in 1:length(sample_names)){
    cell_names <- rownames(metadata[metadata$CellType.Corrected == cell_type_names[j] &
                                      metadata$Sample == sample_names[k], ])
    if (length(cell_names) > 3){
      markers_count[, i] <- rowMeans(markers_flt_count[, cell_names])
      colnames_counts <- c(colnames_counts, paste(cell_type_names[j], sample_names[k], sep = "-"))
      i = i+1
    }
  }
}
markers_count <- as.matrix(markers_count)
markers_count <- markers_count[,colSums(markers_count, na.rm = T) != 0]
markers_count <- markers_count[!is.na(markers_count[,1]),]
colnames(markers_count) <- colnames_counts

markers_count_bis <- flexible_normalization(markers_count+1, by_row = T)
markers_count_bis <- flexible_normalization(markers_count+1, by_row = F)
markers_count[markers_count > 4] <- 4
markers_count[markers_count < -1.5] <- -1.5

annotation <- data.frame(CellType = factor(unlist(strsplit(colnames(markers_count), 
                                              split = "-", fixed = T))[seq(1,2*length(colnames_counts),2)], levels = cell_type_names))
rownames(annotation) <- colnames(markers_count)
ann_colors <- list(CellType = cell_type_color)

pp <- pheatmap::pheatmap(markers_count, cluster_rows = F, cluster_cols = F, 
                         annotation_col = annotation, annotation_colors = ann_colors,     
                         angle_col = 90, annotation_legend = F)

# save_pheatmap_pdf(pp, "Heatmap_markers_samples.pdf", width = 40, height = 40)


#### EPITHELIAL NASAL-------------------------------------------------------------------------

epithelial <- c( 'Suprabasal N', 'Secretory N',  'Multiciliated N')

markers <- c()
for (i in epithelial) {
  markers <- c(markers, unlist(marker_genes[,i]))
}

markers <- markers[!is.na(markers)]
markers <- markers[markers %in% rownames(countTable)]
markers_flt_count <- countTable[markers, ]
# colnames(markers_flt_count) <- gsub(".", "-", colnames(countTable), fixed = T)

cell_type_names <- c("Cycling Basal", "Basal", "Suprabasal", "Suprabasal N", "Secretory", 
                     "Secretory N", "Deuterosomal", "Multiciliated", "Multiciliated N",
                     "SMG Goblet", "Serous",  "AT1", "AT2",  'Brush cells', 'PNEC', 'Ionocyte')
sample_names <- unique(metadata[order(metadata$Position, metadata$Method), "Sample"])
markers_count <- matrix(0, nrow = length(markers), ncol = length(cell_type_names)*length(unique(metadata$Sample)))

rownames(markers_count) <- markers
colnames_counts <- c()
i=1
for (j in 1:length(cell_type_names)){
  print(cell_type_names[j])
  for (k in 1:length(sample_names)){
    cell_names <- rownames(metadata[metadata$CellType.Corrected == cell_type_names[j] &
                                      metadata$Sample == sample_names[k], ])
    if (length(cell_names) > 3){
      markers_count[, i] <- rowMeans(markers_flt_count[, cell_names])
      colnames_counts <- c(colnames_counts, paste(cell_type_names[j], sample_names[k], sep = "-"))
      i = i+1
    }
  }
}
markers_count <- as.matrix(markers_count)
markers_count <- markers_count[,colSums(markers_count, na.rm = T) != 0]
markers_count <- markers_count[!is.na(markers_count[,1]),]
colnames(markers_count) <- colnames_counts

markers_count_bis <- flexible_normalization(markers_count+1, by_row = T)
markers_count_bis <- flexible_normalization(markers_count+1, by_row = F)
markers_count[markers_count > 2] <- 2
markers_count[markers_count < -1.5] <- -1.5

annotation <- data.frame(CellType = factor(unlist(strsplit(colnames(markers_count), 
                                                           split = "-", fixed = T))[seq(1,2*length(colnames_counts),2)], levels = cell_type_names))
rownames(annotation) <- colnames(markers_count)
ann_colors <- list(CellType = cell_type_color)

pp <- pheatmap::pheatmap(markers_count, cluster_rows = F, cluster_cols = F, 
                         annotation_col = annotation, annotation_colors = ann_colors,     
                         angle_col = 90, annotation_legend = F)

# save_pheatmap_pdf(pp, "Heatmap_markers_samples_nasal.pdf", width = 40, height = 40)


#### IMMUNE -------------------------------------------------------------------------

markers <- c()
for (i in immune) {
  markers <- c(markers, unlist(marker_genes[,i]))
}

markers <- markers[!is.na(markers)]
markers <- markers[markers %in% rownames(countTable)]
markers_flt_count <- countTable[markers, ]
# colnames(markers_flt_count) <- gsub(".", "-", colnames(countTable), fixed = T)

cell_type_names <- c( 'B cells', 'Plasma cells', 'LT/NK', 'Mast cells', 'Dendritic', 'Monocyte', 'Macrophage')
sample_names <- unique(metadata[order(metadata$Position, metadata$Method), "Sample"])
markers_count <- matrix(0, nrow = length(markers), ncol = length(cell_type_names)*length(unique(metadata$Sample)))

rownames(markers_count) <- markers
colnames_counts <- c()
i=1
for (j in 1:length(cell_type_names)){
  print(cell_type_names[j])
  for (k in 1:length(sample_names)){
    cell_names <- rownames(metadata[metadata$CellType.Corrected == cell_type_names[j] &
                                      metadata$Sample == sample_names[k], ])
    if (length(cell_names) > 3){
      markers_count[, i] <- rowMeans(markers_flt_count[, cell_names])
      colnames_counts <- c(colnames_counts, paste(cell_type_names[j], sample_names[k], sep = "-"))
      i = i+1
    }
  }
}
markers_count <- as.matrix(markers_count)
markers_count <- markers_count[,colSums(markers_count, na.rm = T) != 0]
markers_count <- markers_count[!is.na(markers_count[,1]),]
colnames(markers_count) <- colnames_counts

markers_count_bis <- flexible_normalization(markers_count+1, by_row = T)
markers_count_bis <- flexible_normalization(markers_count+1, by_row = F)
markers_count[markers_count > 4] <- 4
markers_count[markers_count < -1.5] <- -1.5

annotation <- data.frame(CellType = factor(unlist(strsplit(colnames(markers_count), 
                                              split = "-", fixed = T))[seq(1,2*length(colnames_counts),2)], levels = cell_type_names))
rownames(annotation) <- colnames(markers_count)
ann_colors <- list(CellType = cell_type_color)

pp <- pheatmap::pheatmap(markers_count, cluster_rows = F, cluster_cols = F, 
                         annotation_col = annotation, annotation_colors = ann_colors,     
                         angle_col = 90, annotation_legend = F)

# save_pheatmap_pdf(pp, "Heatmap_markers_immune.pdf", width = 40, height = 40)

#### MESENCHYMAL -------------------------------------------------------------------------

markers <- c()
for (i in stromal) {
  markers <- c(markers, unlist(marker_genes[,i]))
}

markers <- markers[!is.na(markers)]
markers <- markers[markers %in% rownames(countTable)]
markers_flt_count <- countTable[markers, ]
# colnames(markers_flt_count) <- gsub(".", "-", colnames(countTable), fixed = T)

cell_type_names <- c( 'Endothelial', 'Fibroblast', 'Pericyte', 'Smooth muscle')
sample_names <- unique(metadata[order(metadata$Position, metadata$Method), "Sample"])
markers_count <- matrix(0, nrow = length(markers), ncol = length(cell_type_names)*length(unique(metadata$Sample)))

rownames(markers_count) <- markers
colnames_counts <- c()
i=1
for (j in 1:length(cell_type_names)){
  print(cell_type_names[j])
  for (k in 1:length(sample_names)){
    cell_names <- rownames(metadata[metadata$CellType.Corrected == cell_type_names[j] &
                                      metadata$Sample == sample_names[k], ])
    if (length(cell_names) > 3){
      markers_count[, i] <- rowMeans(markers_flt_count[, cell_names])
      colnames_counts <- c(colnames_counts, paste(cell_type_names[j], sample_names[k], sep = "-"))
      i = i+1
    }
  }
}
markers_count <- as.matrix(markers_count)
markers_count <- markers_count[,colSums(markers_count, na.rm = T) != 0]
markers_count <- markers_count[!is.na(markers_count[,1]),]
colnames(markers_count) <- colnames_counts

markers_count_bis <- flexible_normalization(markers_count+1, by_row = T)
markers_count_bis <- flexible_normalization(markers_count+1, by_row = F)
markers_count[markers_count > 5] <- 5
markers_count[markers_count < -1.5] <- -1.5

annotation <- data.frame(CellType = factor(unlist(strsplit(colnames(markers_count), 
                                                           split = "-", fixed = T))[seq(1,2*length(colnames_counts),2)], levels = cell_type_names))
rownames(annotation) <- colnames(markers_count)
ann_colors <- list(CellType = cell_type_color)

pp <- pheatmap::pheatmap(markers_count, cluster_rows = F, cluster_cols = F, 
                         annotation_col = annotation, annotation_colors = ann_colors,     
                         angle_col = 90, annotation_legend = F)

# save_pheatmap_pdf(pp, "Heatmap_markers_stromal.pdf", width = 20, height = 20)







