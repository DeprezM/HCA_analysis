library(ggplot2)
library(tidyr)
library(scales)
library(forcats)
library(matrixStats)
library(fgsea)

setwd("/HCA_analysis/6_Figures/Figure_rare_cells")

##########################################################################################
#                       Load datasets & Create variables                                 #
##########################################################################################

cell_type_order <- c('Cycling Basal',
                     'Basal',
                     'Suprabasal', 'Suprabasal N',
                     'Secretory', 'Secretory N',
                     'Deuterosomal', 
                     'Multiciliated', 'Multiciliated N',
                     'AT1', 'AT2', "Mucous Multiciliated cells",
                     'SMG Goblet', 'Serous', 'Rare cells', 
                     'Brush cells', 'PNEC', 'Precursor', 'Ionocyte',
                     'Endothelial', 'Fibroblast', 'Pericyte', 'Smooth muscle',
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


## Focus on rare cells -----------------------------------------------------------------

meta_rare <- metadata[metadata$CellType.Corrected_rare == "Rare cells",]

rare_coords <- read.table(file = "/Data/Focus_rare_cells_umap_coord.tsv", 
                          sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
meta_rare$x_rareCoord <- rare_coords[rownames(meta_rare), 'X0']
meta_rare$y_rareCoord <- rare_coords[rownames(meta_rare), 'X1']
# AT1 and AT2 cells are low abundant cells in this dataset but they are not so-called 'rare cells'
# meta_rare = meta_rare[!meta_rare$CellType.Corrected %in% c("AT1", "AT2"),]

##########################################################################################
#                                        Plots                                           #
##########################################################################################

### UMAP ---------------------------------------------------------------------------------

ggplot(meta_rare, aes(x = x_rareCoord, y = y_rareCoord, fill = CellType.Corrected))+
  geom_point(shape = 21, colour = "grey20", size = 3, stroke = 0.4) + 
  scale_fill_manual(values = cell_type_color) +
  theme_classic() +
  theme(legend.position = "none")

# pdf("umap_rare_cell_type.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()

### Proportion ---------------------------------------------------------------------------

df = as.data.frame.matrix(table(meta_rare$CellType.Corrected, meta_rare$Position))
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

# pdf("proportion_rare_position.pdf", width=6, height=5, useDingbats=FALSE)
# dev.off()


df = as.data.frame.matrix(table(meta_rare$CellType.Corrected, meta_rare$Method))
#write.csv(df, file = "CellType_composition.csv", quote = F)
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

# pdf("proportion_rare_method.pdf", width=6, height=5, useDingbats=FALSE)
# dev.off()


df = as.data.frame.matrix(table(metadata$CellType.Corrected, metadata$Position))
df = as.data.frame(df / rowSums(df))
df$cell_type <- rownames(df)
df_long = gather(df, position, percentage, 'Nasal':'Distal')
df_long$position <- factor(df_long$position, levels = names(position_color))


ggplot(df_long[df_long$cell_type == "Mucous Multiciliated cells",], aes(x = "", y = percentage, fill = position)) +
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

# pdf("proportion_doublePos_position.pdf", width=3, height=3.5, useDingbats=FALSE)
# dev.off()


df = as.data.frame.matrix(table(metadata$CellType.Corrected, metadata$Method))
df = as.data.frame(df / rowSums(df))
df$cell_type <- rownames(df)
df_long = gather(df, position, percentage, 'Biopsy':'Brushing')
df_long$position <- factor(df_long$position, levels = names(method_color))

ggplot(df_long[df_long$cell_type == "Mucous Multiciliated cells",], aes(x = "", y = percentage, fill = position)) +
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

# pdf("proportion_doublePos_method.pdf", width=3, height=3.5, useDingbats=FALSE)
# dev.off()

### Dot Plot  ------------------------------------------------------------------------

meta_rare <- metadata[metadata$CellType.Corrected %in% c("Mucous Multiciliated cells",
                                                         'Brush cells', 'PNEC', 'Precursor', 'Ionocyte'),]

load('/Data/scranNorm_dataset.Rda')
colnames(countTable) <- gsub(".", "-", colnames(countTable), fixed = T)


markers <- c(#"SFTPA2", "SFTPB", "SFTPC", "SFTPD", 
             #"LRRK2", "CA2", "SLC34A2", "NAPSA", "HOPX",# AT2
             
             #"SFTA2", "AGER", "SLC39A8", "TNNC1", "RTKN2", 
             #"ANXA3", "CAV2", "SPOCK2", "CLDN18", "CAV1", # AT1
  
             "IGFBP2", "ABRACL", "MYB", "NREP", "MYC" , "KYNU",              
             "MARCKSL1", "BIK", "AZGP1",'CRYM',"HES6","MDK","LRMP",  # Precursor
             
             "ASCL3", "RARRES2", "STAP1", "TFCP2L1", "CLCNKB", 
             "ATP6V1G3", "TMEM61", "CLCNKA", "CFTR", "FOXI1", # Ionocyte
             
             "RGS13", "PBXIP1", "BMX", "PLCG2",
             "C15orf48", "AVIL", "HCK", "C11orf53", # Brush cells            

             "GRP", "PCSK1N", "SCG2", "SCGN", "IGFBP5",
             "MS4A8", "NEB", "CPE", "MIAT", "BEX1", # PNEC
             
             "FOXJ1", "PIFO","TPPP3", "MUC5AC", "VMO1", "BPIFA1", "MSMB", "PSCA" #DoublePos
)

colnames(countTable) <- gsub(".", "-", colnames(countTable), fixed = T)
flt_countTable <- countTable[markers, rownames(meta_rare)]

df_meanExpr <- as.data.frame(matrix(0, ncol = 1+length(unique(meta_rare$CellType.Corrected)), 
                                    nrow = length(markers)))
df_perc <- as.data.frame(matrix(0, ncol = 1 + length(unique(meta_rare$CellType.Corrected)), 
                                nrow = length(markers)))
colnames(df_meanExpr) <- colnames(df_perc) <- c("Gene", unique(meta_rare$CellType.Corrected))
df_meanExpr$Gene <- df_perc$Gene <- rownames(flt_countTable)
for (j in unique(meta_rare$CellType.Corrected)){
  cells_id = rownames(meta_rare[meta_rare$CellType.Corrected == j,])
  df_meanExpr[, j] <- rowMeans(flt_countTable[, cells_id])
  df_perc[, j] <- (rowSums(flt_countTable[, cells_id] > 0) / length(cells_id)) *100
}

meanExpr_long <- gather(df_meanExpr, CellType, Expression, "Precursor":"Brush cells")
perc_long <- gather(df_perc, CellType, Percentage, "Precursor":"Brush cells")

dotplot_df <- cbind(meanExpr_long, Percentage = perc_long$Percentage)
# dotplot_df$Gene <- factor(dotplot_df$Gene, levels = markers)
dotplot_df$CellType <- factor(dotplot_df$CellType, 
                              levels = c("Mucous Multiciliated cells","PNEC", "Brush cells",  "Ionocyte", "Precursor"))

dotplot_df <- dotplot_df[!is.na(dotplot_df$Expression),]
dotplot_df$Expression[dotplot_df$Expression > 4] <- 4
ggplot(dotplot_df, aes(y = CellType, x = fct_inorder(Gene), size = rescale(Percentage),
                       alpha = rescale(Expression))) +
  geom_point(color = "#1a3365") +
  theme_classic() +
  scale_size(name = "Percentage \nExpressing") +
  scale_alpha(name = "Mean \nExpression") +#range = c(0.3, 1),
  # scale_color_manual(values = colors_cluster, guide=F ) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        axis.text = element_text(size = 10, face = "bold"))

# pdf("markers_rare_cells_test_2.pdf", width=14, height=4, useDingbats=FALSE)
# dev.off()

### Dot Plot (all epithelial cells)  ---------------------------------------------------------

load('/Data/SoupX_norm_dataset.Rda')
markers = markers[markers %in% rownames(countTable)]

cell_type_to_keep <- c("Basal", "Suprabasal", "Secretory", "Multiciliated", "Deuterosomal", 
                       "Suprabasal N", "Secretory N","Multiciliated N",
                       "Mucous Multiciliated cells",
                       'Brush cells', 'PNEC', 'Precursor', 'Ionocyte')

flt_countTable <- countTable[markers, 
                             rownames(metadata[metadata$CellType.Corrected %in% cell_type_to_keep, ])]

metadata <- metadata[metadata$CellType.Corrected %in% cell_type_to_keep, ]
df_meanExpr <- as.data.frame(matrix(0, ncol = 1+length(unique(metadata$CellType.Corrected)), 
                                    nrow = length(markers)))
df_perc <- as.data.frame(matrix(0, ncol = 1 + length(unique(metadata$CellType.Corrected)), 
                                nrow = length(markers)))
colnames(df_meanExpr) <- colnames(df_perc) <- c("Gene", unique(metadata$CellType.Corrected))
df_meanExpr$Gene <- df_perc$Gene <- rownames(flt_countTable)
for (j in unique(metadata$CellType.Corrected)){
  cells_id = rownames(metadata[metadata$CellType.Corrected == j,])
  df_meanExpr[, j] <- rowMeans(flt_countTable[, cells_id])
  df_perc[, j] <- (rowSums(flt_countTable[, cells_id] > 0) / length(cells_id)) *100
}

meanExpr_long <- gather(df_meanExpr, CellType, Expression, "Suprabasal N":"Brush cells")
perc_long <- gather(df_perc, CellType, Percentage, "Suprabasal N":"Brush cells")

dotplot_df <- cbind(meanExpr_long, Percentage = perc_long$Percentage)
dotplot_df$CellType <- factor(dotplot_df$CellType, 
                              levels = c("Basal", "Suprabasal", "Suprabasal N", "Secretory", "Secretory N",
                                         "Deuterosomal", "Multiciliated", "Multiciliated N", 
                                         'Precursor', 'Ionocyte', 'Brush cells', 'PNEC', "Mucous Multiciliated cells"))

dotplot_df <- dotplot_df[!is.na(dotplot_df$Expression),]
dotplot_df$Expression[dotplot_df$Expression > 4] <- 4
ggplot(dotplot_df, aes(y = CellType, x = fct_inorder(Gene), size = rescale(Percentage),
                       alpha = rescale(Expression))) +
  geom_point(color = "#1a3365") +
  
  theme_classic() +
  scale_size(name = "Percentage \nExpressing") +
  scale_alpha(name = "Mean \nExpression") +#range = c(0.3, 1),
  # scale_color_manual(values = colors_cluster, guide=F ) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        axis.text = element_text(size = 10, face = "bold"))

# pdf("markers_rare_cells_all_ct_test.pdf", width=20, height=8, useDingbats=FALSE)
# dev.off()


### Dot Plot - AT1 AT2 --------------------------------------------------------------------

markers <- c("SFTPA2", "SFTPB", "SFTPC", "SFTPD", 
  "LRRK2", "CA2", "SLC34A2", "NAPSA", "HOPX",# AT2
  
  "SFTA2", "AGER", "SLC39A8", "TNNC1", "RTKN2", 
  "ANXA3", "CAV2", "SPOCK2", "CLDN18", "CAV1" # AT1
  )

markers = markers[markers %in% rownames(countTable)]
cell_type_to_keep <- c("Basal", "Suprabasal", "Secretory", "Multiciliated", "Deuterosomal", 
                       "Suprabasal N", "Secretory N","Multiciliated N",
                       "AT1",
                       'AT2')

flt_countTable <- countTable[markers, 
                             rownames(metadata[metadata$CellType.Corrected %in% cell_type_to_keep, ])]

metadata <- metadata[metadata$CellType.Corrected %in% cell_type_to_keep, ]
df_meanExpr <- as.data.frame(matrix(0, ncol = 1+length(unique(metadata$CellType.Corrected)), 
                                    nrow = length(markers)))
df_perc <- as.data.frame(matrix(0, ncol = 1 + length(unique(metadata$CellType.Corrected)), 
                                nrow = length(markers)))
colnames(df_meanExpr) <- colnames(df_perc) <- c("Gene", unique(metadata$CellType.Corrected))
df_meanExpr$Gene <- df_perc$Gene <- rownames(flt_countTable)
for (j in unique(metadata$CellType.Corrected)){
  cells_id = rownames(metadata[metadata$CellType.Corrected == j,])
  df_meanExpr[, j] <- rowMeans(flt_countTable[, cells_id])
  df_perc[, j] <- (rowSums(flt_countTable[, cells_id] > 0) / length(cells_id)) *100
}

meanExpr_long <- gather(df_meanExpr, CellType, Expression, "Suprabasal N":"AT1")
perc_long <- gather(df_perc, CellType, Percentage, "Suprabasal N":"AT1")

dotplot_df <- cbind(meanExpr_long, Percentage = perc_long$Percentage)
# dotplot_df$Gene <- factor(dotplot_df$Gene, levels = markers)
dotplot_df$CellType <- factor(dotplot_df$CellType, 
                              levels = c("Basal", "Suprabasal", "Suprabasal N", "Secretory", "Secretory N",
                                         "Deuterosomal", "Multiciliated", "Multiciliated N", 
                                         "AT1", "AT2"))

dotplot_df <- dotplot_df[!is.na(dotplot_df$Expression),]
dotplot_df$Expression[dotplot_df$Expression > 4] <- 4
ggplot(dotplot_df, aes(y = CellType, x = fct_inorder(Gene), size = rescale(Percentage),
                       alpha = rescale(Expression))) +
  geom_point(color = "#1a3365") +
  
  theme_classic() +
  scale_size(name = "Percentage \nExpressing") +
  scale_alpha(name = "Mean \nExpression") +#range = c(0.3, 1),
  # scale_color_manual(values = colors_cluster, guide=F ) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        axis.text = element_text(size = 10, face = "bold"))

# pdf("markers_AT1_AT2.pdf", width=6, height=5, useDingbats=FALSE)
# dev.off()



### Heatmap TF -----------------------------------------------------------------------------

load("/Data/TFscore_matrix.Rda")
# load("/home/truchi/HumanCellAtlas/TF_score/TFscore_TopTable.Rda")
# write.table(TF_matrix_Seurat.markers, 
#             file = '/data/deprez_data/HCA/Analysis/FullDataset_v4/TF_marker.tsv', 
#             sep= "\t", quote = F)
# 

tf_matrix <- as.matrix(TF_matrix@assays$RNA@counts)

markers <- c(#"RTKN2", "AGER", "RTKN2", "IRX2",  # AT1
             
             #"ETV5", "MLXIPL", "SOX5", "ZNF77", "ZFP28",  # AT2 
  
             "POU2F3", "SOX4", "MYC", "ASCL2", # Precursor
             
             "ASCL3", "FOXI1", "DMRT2","SRF", "BNC2", # Ionocyte
             
             "HOXC5", "HMX2", "ANXA4","ARNT2", "OVOL3",   # Brush cells
            
             "NR0B1","HOXB5", "ASCL1", "INSM1","FOXA2" # PNEC
)

markers_flt_count <- tf_matrix[markers, ]
colnames(markers_flt_count) <- gsub(".", "-", colnames(countTable), fixed = T)

markers_count <- matrix(0, nrow = length(markers), ncol = 4)#length(cell_type_order) - 2
colnames(markers_count) <- c("Precursor","Ionocyte", "Brush cells",
                             "PNEC")

rownames(markers_count) <- markers

for (i in 1:ncol(markers_count)){
  print(colnames(markers_count)[i])
  cell_names <- rownames(metadata[metadata$CellType.Corrected == colnames(markers_count)[i], ])
  markers_count[, colnames(markers_count)[i]] <- rowMedians(markers_flt_count[, cell_names])
}
markers_count <- log(markers_count+ 1)

annotation <- data.frame(CellType = factor(colnames(markers_count), levels = cell_type_order))
rownames(annotation) <- colnames(markers_count)
ann_colors <- list(CellType = cell_type_color)

pp <- pheatmap::pheatmap(markers_count, cluster_rows = F, cluster_cols = F, 
                         annotation_col = annotation, annotation_colors = ann_colors,     
                         angle_col = 90, annotation_legend = F)

# save_pheatmap_pdf(pp, "Heatmap_markers_TF_rare_cells.pdf", width = 5, height = 6)



### Correlation Heatmap ---------------------------------------------------------------------

load('/Data/scranNorm_dataset.Rda')
colnames(countTable) <- gsub(".", "-", colnames(countTable), fixed = T)

meta_corr <- metadata
meta_corr$nasal_celltype <- as.character(metadata$CellType.Corrected)

celltype <- c("Basal", "Suprabasal", "Secretory", "Multiciliated", "Deuterosomal", 
              "Suprabasal N", "Secretory N", "Multiciliated N", 
              "Brush cells", "Ionocyte", "PNEC", "Precursor")
meta <- meta_corr[meta_corr$nasal_celltype %in% celltype, ]

celltype_B <- c("Basal", "Suprabasal", "Secretory", "Multiciliated", "Deuterosomal")
meta$position_cell_type <- meta$nasal_celltype
meta$position_cell_type[(meta$nasal_celltype %in% celltype_B) & (meta$Position == "Nasal")] <- 
  paste0(meta$position_cell_type[(meta$nasal_celltype %in% celltype_B) & (meta$Position == "Nasal")], " N")
meta$position_cell_type[(meta$nasal_celltype %in% celltype_B) & (meta$Position == "Proximal")] <- 
  paste0(meta$position_cell_type[(meta$nasal_celltype %in% celltype_B) & (meta$Position == "Proximal")], " P")
meta$position_cell_type[(meta$nasal_celltype %in% celltype_B) & (meta$Position == "Intermediate")] <- 
  paste0(meta$position_cell_type[(meta$nasal_celltype %in% celltype_B) & (meta$Position == "Intermediate")], " I")
meta$position_cell_type[(meta$nasal_celltype %in% celltype_B) & (meta$Position == "Distal")] <- 
  paste0(meta$position_cell_type[(meta$nasal_celltype %in% celltype_B) & (meta$Position == "Distal")], " D")


table(meta$position_cell_type)
celltype = unique(meta$position_cell_type)

norm = countTable[, rownames(meta)]
dim(norm)

cmed = c()
for (i in 1:length(celltype)) {
  print(celltype[i])
  cmed = cbind(cmed, rowMeans(norm[, meta$position_cell_type==celltype[i]]))
}
colnames(cmed)=celltype
ct = cor(cmed)
unlist_el <- unlist(strsplit(colnames(ct), " "))


cell_type_column <- c()
position_column <- c()
for (i in celltype) {
  cell_type_column <- c(cell_type_column, unlist(strsplit(i, " "))[1])
  if (length(unlist(strsplit(i, " "))) == 2){
    position_column <- c(position_column, unlist(strsplit(i, " "))[2])
  } else {
    position_column <- c(position_column, " ")
  }
}

annotation <- data.frame(cell_type = cell_type_column,
                         position = position_column,
                         stringsAsFactors = F)
rownames(annotation) <- colnames(ct)
annotation["Brush cells", "position"] <- " "
annotation["Brush cells", "cell_type"] <- "Brush cells"
annotation["Secretory N", "cell_type"] <- "Secretory N"
annotation["Suprabasal N", "cell_type"] <- "Suprabasal N"
annotation["Multiciliated N", "cell_type"] <- "Multiciliated N"

position_color <- c('N' = '#ffd966',
                    'P' = '#e69138',
                    'I' = '#e06666',
                    'D' = '#85200c',
                    ' ' = '#b3b3b3')


anno_col <- list(cell_type = cell_type_color, position = position_color)

pheatmap::pheatmap(ct, clustering_method = "ward.D")
pp <- pheatmap::pheatmap(ct,  main="", 
               annotation_col = annotation, annotation_colors = anno_col)

# save_pheatmap_pdf(pp, "Heatmap_correlation_epithelial_rare_cell_type.pdf", width = 8.5, height = 7)



### Umap Multiciliating Goblet -------------------------------------------------------------

load('/Data/SoupX_norm_dataset.Rda')

metadata_genes <- metadata
metadata_genes$Gene1 <- countTable["MUC5AC",rownames(metadata_genes)]
metadata_genes$Gene2 <- countTable["FOXJ1", rownames(metadata_genes)]
metadata_genes$coExpr <- "none"
metadata_genes$coExpr[metadata_genes$Gene1 > 0] <- "Gene 1"
metadata_genes$coExpr[metadata_genes$Gene2 > 0] <- "Gene 2"
metadata_genes$coExpr[metadata_genes$Gene1 > 0 & metadata_genes$Gene2 > 0] <- "Coexpr"
metadata_genes$isCo <- 1
metadata_genes$isCo[metadata_genes$coExpr == "Coexpr"] <- 2

metadata_genes$coExpr <- factor(metadata_genes$coExpr, 
                                levels = c("none", "Gene 1", "Gene 2", "Coexpr"))

ggplot(metadata_genes, 
       aes(x = x, y = y, color = coExpr, size = isCo, shape = (isCo == 2))) +
  geom_point() +
  geom_point(aes(x = x, y = y, color = coExpr, size = isCo, shape = (isCo == 2)), 
             data = subset(metadata_genes, coExpr == 'Coexpr')) +
  scale_color_manual(values = c("Gene 1" = "#006600", "Gene 2" = "#466cb9", "Coexpr" = "#ff00ff", "none" = "grey")) +
  scale_shape_manual(values = c(19, 17)) +
  scale_size(range = c(1.5, 3)) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill = "transparent"))


table(metadata_genes$coExpr, metadata_genes$Position)

# pdf("umap_double_pos_SOUPX.pdf", width=4, height=4, useDingbats=FALSE)
# dev.off()

