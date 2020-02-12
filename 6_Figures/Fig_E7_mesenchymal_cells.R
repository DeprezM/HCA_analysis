library(ggplot2)
library(tidyr)
library(scales)
library(forcats)
library(matrixStats)


setwd("/HCA_analysis/6_Figures/Figure_mesenchymal_cells")

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
                     'Endothelial', 'Fibroblast', 'Myofibroblast', 'Smooth muscle',
                     'B cells', 'Mast cells', 'LT/NK', 'Plasma cells','Monocyte', 'Dendritic', 'Macrophage')

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

## Focus on mesenchymal cells -----------------------------------------------------------------
resident_cell_type <- c("Endothelial", "Smooth muscle", "Pericyte", "Fibroblast")
meta_resident <- metadata[metadata$CellType.Corrected %in% resident_cell_type,]

resident_coord <- read.table(file = "/Data/Focus_mesenchymal_cells_metadata.tsv",
                                     sep = "\t", header = T, row.names = 1)
meta_resident$x_rareCoord <- resident_coord[rownames(meta_resident), 'umap_1']
meta_resident$y_rareCoord <- resident_coord[rownames(meta_resident), 'umap_2']



##########################################################################################
#                                        Plots                                           #
##########################################################################################

### UMAP ---------------------------------------------------------------------------------

ggplot(meta_resident, aes(x = x_rareCoord, y = y_rareCoord, color = CellType.Corrected))+
  geom_point(size = 1) + 
  scale_color_manual(values = cell_type_color) +
  theme_classic() +
  theme(legend.position = "none")

# pdf("umap_resident_cell_type.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()

### Proportion (Pie charts) ---------------------------------------------------------------

df = as.data.frame.matrix(table(meta_resident$CellType.Corrected, meta_resident$Position))
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

# pdf("proportion_resident_position.pdf", width=6, height=6, useDingbats=FALSE)
# dev.off()


df = as.data.frame.matrix(table(meta_resident$CellType.Corrected, meta_resident$Method))
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
# 
# pdf("proportion_rare_method.pdf", width=6, height=6, useDingbats=FALSE)
# dev.off()


### Proportion (bar chart) ---------------------------------------------------------------
df = as.data.frame.matrix(table(metadata$Sample, metadata$CellType.Corrected_rare))
df = as.data.frame(df / rowSums(df))
df$Sample <- factor(rownames(df), 
                    levels = unique(metadata[order(metadata$SampleType, metadata$Position, metadata$Method), "Sample"]))
df_long = gather(df, cell_type, percentage, 'B cells':'Suprabasal N')
df_long$cell_type <- factor(df_long$cell_type, levels = cell_type_order)


ggplot(df_long[df_long$cell_type %in% c('Endothelial', 'Fibroblast', 'Myofibroblast', 'Smooth muscle'),], 
       aes(x = Sample, y = percentage, fill = cell_type)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = cell_type_color) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, angle = 90),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size =14),
    axis.text.y = element_text(size = 14),
    axis.line = element_line(colour = "grey50"),
    axis.ticks = element_line(colour = "grey50"),
    legend.position = "none" 
  )

# pdf("proportion_cell_type_per_sample_cycling_mesenchyme_1.pdf", width = 8, height = 3, useDingbats=FALSE)
# dev.off()



### Heatmap TF activity -------------------------------------------------------------------
load("/Data/TFscore_matrix.Rda")
# load("/home/truchi/HumanCellAtlas/TF_score/TFscore_TopTable.Rda")
# write.table(TF_matrix_Seurat.markers, 
#             file = '/data/deprez_data/HCA/Analysis/FullDataset_v4/TF_marker.tsv', 
#             sep= "\t", quote = F)


tf_matrix <- as.matrix(TF_matrix@assays$RNA@counts)

markers <- c("CD34", "MEOX1", "LDB2","SOX18", "ERG", # Endothelial
             
             "OSR2", "NR2F1", "GLI1", "SHOX2", "FOXF1", # Fibroblast
             
             "DUXA", "HIF3A", "HOXA4", "ABRA", "FLNA",  # Smooth muscle
             
             "EBF2",
             "EBF1",
             "HOXA7",
             "ZNF660","HEY2"# Myofibroblast
)

markers %in% rownames(tf_matrix)
markers_flt_count <- tf_matrix[markers, ]
colnames(markers_flt_count) <- gsub(".", "-", colnames(countTable), fixed = T)

markers_count <- matrix(0, nrow = length(markers), ncol = 4)#length(cell_type_order) - 2
colnames(markers_count) <- c("Endothelial", "Fibroblast",
                             "Smooth muscle", "Myofibroblast")

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

# save_pheatmap_pdf(pp, "Heatmap_markers_TF_resident.pdf", width = 5, height = 6)


### Violin plot for Gene expression ---------------------------------------------------

load('/Data/scranNorm_dataset.Rda')
colnames(countTable) <- gsub(".", "-", colnames(countTable), fixed = T)


gene <- c("MYH11", "PPP1R14A", "MYLK",  "MYH10",  "ACTN1", "CNN1",
          "PDLIM3", "CKB", "RAMP1",
          "CRIP1", "RERGL", "IGFBP5", "NDUFA4L2", "NOTCH3", "FRZB", "RGS5","MDK")

metadata_genes <- cbind(metadata, t(countTable[gene,]))
meta_subset <- metadata_genes[metadata_genes$CellType.Corrected %in% c("Smooth muscle", "Myofibroblast", "Fibroblast"),]
meta_subset <- meta_subset[, c("CellType.Corrected", gene)]
df_long <- gather(meta_subset, Gene, Value, gene)
df_long$Gene <- factor(df_long$Gene, levels = gene)
head(df_long)

ggplot(df_long, 
       aes(x = CellType.Corrected, y = Value, fill = CellType.Corrected)) +
  geom_violin(scale = "width") + 
  scale_fill_manual(values = cell_type_color)+
  scale_y_continuous(position = "right", breaks = c(min(df_long$Value), max(df_long$Value)))+
  facet_wrap(~Gene, ncol = 1, strip.position="left", scales = "free_y") +
  theme_classic()+
  theme(
    strip.background = element_rect(colour="white", fill="white"),
    axis.text.x = element_text(angle = 90), 
    strip.text.y = element_text(angle = 180),
    # axis.text.y = element_text(size = 15),
    axis.title = element_blank())

# 
# pdf(paste0("Violin_smooth musccle.pdf"), width = 5, height = 9, useDingbats=FALSE)
# dev.off()
