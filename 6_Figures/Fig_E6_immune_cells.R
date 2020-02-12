library(ggplot2)
library(tidyr)
library(scales)
library(forcats)
library(matrixStats)

setwd("/HCA_analysis/6_Figures/Figure_immune_cells")

##########################################################################################
#                       Load datasets & Create variables                                 #
##########################################################################################

## Main variables and functions ----------------------------------------------------------
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

sample_type_color <- c('Nasal\nBrushing' = '#ffd966',
                       'Nasal\nBiopsy' = '#cc9900',
                       'Proximal\nBiopsy' = '#e69138',
                       'Intermediate\nBiopsy' = '#e06666',
                       'Distal\nBrushing' = '#85200c')

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
colnames(metadata)

## Focus on immune cells -----------------------------------------------------------------
immune_cell_type <- c("Dendritic", "B cells", "Plasma cells", "Mast cells", 
                      "Macrophage", "LT/NK", "Monocyte")
meta_immune <- metadata[metadata$CellType.Corrected %in% immune_cell_type,]

immune_coord <- read.table(file = "/Data/Focus_immune_cells_metadata.tsv", 
                           sep = "\t", header = T, stringsAsFactors = F, row.names = 1)

meta_immune$x_rareCoord <- immune_coord[rownames(meta_immune), 'umap_1']
meta_immune$y_rareCoord <- immune_coord[rownames(meta_immune), 'umap_2']


##########################################################################################
#                                        Plots                                           #
##########################################################################################

### UMAP ---------------------------------------------------------------------------------

ggplot(meta_immune, aes(x = x_rareCoord, y = y_rareCoord, color = CellType.Corrected))+
  geom_point(size = 0.5) + 
  scale_color_manual(values = cell_type_color) +
  theme_classic() #+
  theme(legend.position = "none")

# pdf("umap_immune_cell_type.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()

### UMAP & Gene expression ---------------------------------------------------------------  
  
load('/Data/scranNorm_dataset.Rda')
colnames(countTable) <- gsub(".", "-", colnames(countTable), fixed = T)

gene = "CD2"
meta_immune$gene <- t(countTable[gene, rownames(meta_immune)])

ggplot(subset(meta_immune, CellType.Corrected == "LT/NK" & gene < 1 & x_rareCoord <0), 
       aes(x = x_rareCoord, y = y_rareCoord))+
  geom_point(size= 1, color = "grey70") +
  geom_point(data = subset(meta_immune, CellType.Corrected == "LT/NK" & gene > 1 & x_rareCoord <0), 
             aes(x = x_rareCoord, y = y_rareCoord, color = gene), size = 1.5) +
  scale_color_gradient(low = "grey70", high = "red4") +
  theme_classic() +
  theme(legend.position = "none")

# pdf(paste0("umap_immune_",gene,".pdf"), width=4, height=2, useDingbats=FALSE)
# dev.off()


### Proportion (Pie charts) ---------------------------------------------------------------

df = as.data.frame.matrix(table(meta_immune$CellType.Corrected, meta_immune$Position))
df = as.data.frame(df / rowSums(df))
df$cell_type <- rownames(df)
df_long = gather(df, position, percentage, 'Nasal':'Distal')
df_long$position <- factor(df_long$position, levels = names(position_color))


ggplot(df_long, aes(x = "", y = percentage, fill = position)) +
  geom_bar(stat = 'identity', color = "grey20") +
  scale_fill_manual(values = position_color) +
  theme_classic() +
  coord_polar("y", start = 0, direction = -1) +
  facet_wrap(~cell_type) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none" 
  )

# pdf("proportion_immune_position.pdf", width=6, height=6, useDingbats=FALSE)
# dev.off()


df = as.data.frame.matrix(table(meta_immune$CellType.Corrected, meta_immune$Method))
df = as.data.frame(df / rowSums(df))
df$cell_type <- rownames(df)
df_long = gather(df, position, percentage, 'Biopsy':'Brushing')
df_long$position <- factor(df_long$position, levels = names(method_color))

ggplot(df_long, aes(x = "", y = percentage, fill = position)) +
  geom_bar(stat = 'identity', color = "grey20") +
  scale_fill_manual(values = method_color) +
  theme_classic() +
  coord_polar("y", start = 0, direction = -1) +
  facet_wrap(~cell_type) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none" 
  )

# pdf("proportion_rare_method.pdf", width=6, height=6, useDingbats=FALSE)
# dev.off()

### Proportion (bar chart) ---------------------------------------------------------------
df = as.data.frame.matrix(table(metadata$Sample, metadata$CellType.Corrected_rare))
df = as.data.frame(df / rowSums(df))
df$Sample <- factor(rownames(df), 
                    levels = unique(metadata[order(metadata$SampleType, metadata$Position, metadata$Method), "Sample"]))
df_long = gather(df, cell_type, percentage, 'B cells':'Suprabasal N')
df_long$cell_type <- factor(df_long$cell_type, levels = cell_type_order)


ggplot(df_long[df_long$cell_type %in% c('B cells', 'Plasma cells', 'LT/NK', 'Mast cells', 'Dendritic', 'Monocyte', 'Macrophage'),], 
       aes(x = Sample, y = percentage, fill = cell_type)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = cell_type_color) +
  theme_classic() +
  ylim(0, 0.3)+
  theme(
    axis.text.x = element_text(size = 12, angle = 90),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size =14),
    axis.text.y = element_text(size = 14),
    axis.line = element_line(colour = "grey50"),
    axis.ticks = element_line(colour = "grey50"),
    legend.position = "none" 
  )

pdf("proportion_cell_type_per_sample_cycling_immune.pdf", width = 8, height = 3, useDingbats=FALSE)
dev.off()



### Heatmap TF activity -------------------------------------------------------------------
load("/Data/TFscore_matrix.Rda")
# load("/home/truchi/HumanCellAtlas/TF_score/TFscore_TopTable.Rda")
# write.table(TF_matrix_Seurat.markers, 
#             file = '/data/deprez_data/HCA/Analysis/FullDataset_v4/TF_marker.tsv', 
#             sep= "\t", quote = F)


tf_matrix <- as.matrix(TF_matrix@assays$RNA@counts)

markers <- c("NR1H3", "PPARG", "CARD16", "NLRC4", "NUPR1", # Macrophage
             
             "MAFB", "ZNF331", "NR4A3", "CREM", "NFKB1", # Monocyte
             
             "HIC1", "FIGLA", "CX3CR1", "HCLS1", "ADORA3", # Dendritic
             
             "BCL11B", "STAT4", "XCL1", "CD3D", "BCL11B", # LT/NK
             
             "GATA1", "TAL1", "NTRK1", "GATA2", "SMYD3", # Mast cells
             
             "MIXL1", "ISG20", "E2F5", "IRF4", "ZNF215", # B cells "ISL2"
             
             "SP140", "SPIB", "POU2F2", "PAX5" # Plasma cells
)

markers %in% rownames(tf_matrix)
markers_flt_count <- tf_matrix[markers, ]
colnames(markers_flt_count) <- gsub(".", "-", colnames(countTable), fixed = T)

markers_count <- matrix(0, nrow = length(markers), ncol = 7)#length(cell_type_order) - 2
colnames(markers_count) <- c("Macrophage", "Monocyte", "Dendritic", "LT/NK", 
                             "Mast cells", "B cells", "Plasma cells")

rownames(markers_count) <- markers

for (i in 1:ncol(markers_count)){
  print(colnames(markers_count)[i])
  cell_names <- rownames(metadata[metadata$CellType.Corrected == colnames(markers_count)[i], ])
  markers_count[, colnames(markers_count)[i]] <- rowMedians(markers_flt_count[, cell_names])
}

markers_count[markers_count > 3] <- 3

annotation <- data.frame(CellType = factor(colnames(markers_count), levels = cell_type_order))
rownames(annotation) <- colnames(markers_count)
ann_colors <- list(CellType = cell_type_color)

pp <- pheatmap::pheatmap(markers_count, cluster_rows = F, cluster_cols = F, 
                         annotation_col = annotation, annotation_colors = ann_colors,     
                         angle_col = 90, annotation_legend = F)

save_pheatmap_pdf(pp, "Heatmap_markers_TF_immune.pdf", width = 5, height = 6)





