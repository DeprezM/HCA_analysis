library(ggplot2)
library(tidyr)
library(scales)
library(forcats)

setwd("/HCA_analysis/6_Figures/Figure_cycling_cells")

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
                     'Endothelial', 'Fibroblast', 'Myofibroblast', 'Smooth muscle',
                     'B cells', 'Mast cells', 'LT/NK', 'Plasma cells','Monocyte', 'Dendritic', 'Macrophage')

method_color <-  c('Brushing' = '#009933',
                   'Biopsy' = '#336699')

position_color <- c('Nasal' = '#ffd966',
                    'Proximal' = '#e69138',
                    'Intermediate' = '#e06666',
                    'Distal' = '#85200c')

sampletype_order <- c('Nasal\nBrushing', 'Nasal\nBiopsy', 
                      'Proximal\nBiopsy', 'Intermediate\nBiopsy',
                      'Distal\nBrushing')

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

## Focus on cycling cells -----------------------------------------------------------------
metadata_cycling <- metadata[metadata$CellType.Corrected == "Cycling Basal",]

coord <- read.table(file = "/Data/Focus_cycling_cells_umap_coord.tsv", 
                    sep = "\t", header = T, row.names = 1)
metadata_cycling$Umap1 <- coord[rownames(metadata_cycling),]$X0
metadata_cycling$Umap2 = coord[rownames(metadata_cycling),]$X1

meta_cc <- read.table(file = "/Data/Focus_cycling_cells_metadata.tsv", 
                      sep = "\t", header = T, row.names = 1)
metadata_cycling$S_score <- meta_cc[rownames(metadata_cycling),]$S_score
metadata_cycling$G2M_score <- meta_cc[rownames(metadata_cycling),]$G2M_score
metadata_cycling$G1_score <- meta_cc[rownames(metadata_cycling),]$G1_score

##########################################################################################
#                                        Plots                                           #
##########################################################################################

### UMAP & Cell Cycle score --------------------------------------------------------------

ggplot(metadata_cycling, aes(x = Umap1, y = Umap2, fill = S_score))+
  geom_point(shape = 21, colour = "grey20", size = 3, stroke = 0.4) + 
  # scale_fill_gradient(low = "grey70", high = "red4") +
  scale_fill_gradient2(low = "blue4", high = "red4") +
  theme_classic() + theme(legend.position = "none")

# pdf("umap_cycling_basal_S_score.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()

ggplot(metadata_cycling, aes(x = Umap1, y = Umap2, fill = G2M_score))+
  geom_point(shape = 21, colour = "grey20", size = 3, stroke = 0.4) + 
  # scale_fill_gradient(low = "grey70", high = "red4") +
  scale_fill_gradient2(low = "blue4", high = "red4") +
  theme_classic() + theme(legend.position = "none")

# pdf("umap_cycling_basal_G2M_score.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()

ggplot(metadata_cycling, aes(x = Umap1, y = Umap2, fill = G1_score + 5))+ 
  geom_point(shape = 21, colour = "grey20", size = 3, stroke = 0.4) + 
  # scale_fill_gradient(low = "grey70", high = "red4") +
  scale_fill_gradient2(low = "blue4", high = "red4") +
  theme_classic() + theme(legend.position = "none")
# 
# pdf("umap_cycling_basal_G1_score.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()



### UMAP & Complete dataset ------------------------------------------------------------

ggplot(subset(metadata, CellType.Corrected != "Cycling Basal"), 
       aes(x = x, y = y))+
  geom_point(size= 0.2, color = "grey70") +
  geom_point(data = subset(metadata, CellType.Corrected == "Cycling Basal"), 
             aes(x = x, y = y), color = "#F3766E", size = 1) +
  # scale_color_manual(values = cell_type_grey) +
  theme_classic() +
  theme(legend.position = "none")

# pdf("umap_cycling_basal_total.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()


uncorrected_coord <- read.table(file = "/Data/Umap_uncorrected_embedding.tsv",
                                sep = "\t", header = T, row.names = 1)
colnames(uncorrected_coord)
metadata$unc_umap1 <- uncorrected_coord$umap_1
metadata$unc_umap2 <- uncorrected_coord$umap_2

ggplot(subset(metadata, CellType.Corrected != "Cycling Basal"), 
       aes(x = unc_umap1, y = unc_umap2))+
  geom_point(size= 0.2, color = "grey70") +
  geom_point(data = subset(metadata, CellType.Corrected == "Cycling Basal"), 
             aes(x = unc_umap1, y = unc_umap2), color = "#F3766E", size = 1) +
  # scale_color_manual(values = cell_type_grey) +
  theme_classic() +
  theme(legend.position = "none")

# pdf("umap_uncorrected_cycling_basal_total.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()


### UMAP & Gene Expression --------------------------------------------------------------

load('/Data/scranNorm_dataset.Rda')
colnames(countTable) <- gsub(".", "-", colnames(countTable), fixed = T)

metadata$KRT13 <- t(countTable["KRT13",])
metadata$KRT4 <- t(countTable["KRT4",])

ggplot(subset(metadata, KRT13 < 1), 
       aes(x = x, y = y))+
  geom_point(size= 0.2, color = "grey70") +
  geom_point(data = subset(metadata, KRT13 > 1), 
             aes(x = x, y = y, color = KRT13), size = 0.5) +
  scale_color_gradient(low = "grey70", high = "red4") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(subset(metadata, KRT4 < 1), 
       aes(x = x, y = y))+
  geom_point(size= 0.2, color = "grey70") +
  geom_point(data = subset(metadata, KRT4 > 1), 
             aes(x = x, y = y, color = KRT4), size = 0.5) +
  scale_color_gradient(low = "grey70", high = "red4") +
  theme_classic() +
  theme(legend.position = "none")


# pdf("umap_krt13.pdf", width=3, height=3, useDingbats=FALSE)
# pdf("umap_krt4.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()


### Violin & Cell Cycle Score --------------------------------------------------------------

metadata$S_score <- meta_cc[rownames(metadata),]$S_score
metadata$G2M_score <- meta_cc[rownames(metadata),]$G2M_score
metadata$G1_score <- meta_cc[rownames(metadata),]$G1_score
metadata$CellType.Corrected_rare <- factor(metadata$CellType.Corrected_rare, levels = cell_type_order)

ggplot(metadata, aes(x = CellType.Corrected_rare, y = S_score, fill = CellType.Corrected_rare))+
  geom_violin(scale = "width") + #geom_jitter(size = 0.1, alpha =0.5) + 
  scale_fill_manual(values = cell_type_color) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))

ggplot(metadata, aes(x = CellType.Corrected_rare, y = G2M_score, fill = CellType.Corrected_rare))+
  geom_violin(scale = "width") + #geom_jitter(size = 0.1, alpha =0.5) + 
  scale_fill_manual(values = cell_type_color) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))


# pdf("Cell_cycle_S_score_all_cell_type.pdf", width=6, height=2, useDingbats=FALSE)
# pdf("Cell_cycle_G2M_score_all_cell_type.pdf", width=6, height=2, useDingbats=FALSE)
# 
# dev.off()

### Violin & Gene Expression --------------------------------------------------------------

metadata$KRT13 <- t(countTable["KRT13",])
ggplot(metadata, aes(x = CellType.Corrected_rare, y = KRT13, fill = CellType.Corrected_rare))+
  geom_violin(scale = "width") + geom_jitter(size = 0.1, alpha =0.5) + 
  scale_fill_manual(values = cell_type_color) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))

metadata$KRT13 <- t(countTable["KRT13",])

ggplot(metadata[metadata$CellType.Corrected %in% c("Suprabasal", "Suprabasal N"), ], aes(x = Sample, y = KRT13))+
  geom_violin(scale = "width", fill = cell_type_color["Suprabasal"]) +  
  scale_x_discrete(limits = unique(metadata[order(metadata$Position, metadata$Method), "Sample"])) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))

# pdf("KRT13_suprabasal.pdf", width=6, height=2, useDingbats=FALSE)
# dev.off()

### Proportion (Pie charts) --------------------------------------------------------------

table(meta_cc$position)
meta_cc$position <- factor(meta_cc$position, levels = names(position_color))
df = as.data.frame(table(meta_cc$position))
df$Freq <- df$Freq / sum(df$Freq)

ggplot(df, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = 'identity', color = "grey20") +
  scale_fill_manual(values = position_color) +
  theme_classic() +
  coord_polar("y", start = 0, direction = -1) +
  # facet_wrap(~cell_type) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none" 
  )

# pdf("proportion_cycling_basal_position.pdf", width=2, height=2, useDingbats=FALSE)
# dev.off()

df = as.data.frame(table( meta_cc$method))
df$Freq <- df$Freq / sum(df$Freq)

ggplot(df, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = 'identity', color = "grey20") +
  scale_fill_manual(values = method_color) +
  theme_classic() +
  coord_polar("y", start = 0, direction = -1) +
  # facet_wrap(~cell_type) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none" 
  )

# pdf("proportion_cycling_basal_method.pdf", width=2, height=2, useDingbats=FALSE)
# dev.off()

### Proportion (Barplot) --------------------------------------------------------------

df = as.data.frame.matrix(table(metadata$Sample, metadata$CellType.Corrected_cycling))
df = as.data.frame(df / rowSums(df))
df$Sample <- factor(rownames(df), 
                    levels = unique(metadata[order(metadata$SampleType, metadata$Position, metadata$Method), "Sample"]))
df_long = gather(df, cell_type, percentage, 'Cycling Basal':'Macrophage')
df_long$cell_type <- factor(df_long$cell_type, levels = cell_type_order)

ggplot(df_long, aes(x = Sample, y = percentage, fill = cell_type)) +
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

# pdf("proportion_cell_type_per_sample_cycling.pdf", width = 8, height = 4, useDingbats=FALSE)
# dev.off()

### Dot Plot & Gene Expression --------------------------------------------------------------

load('/Data/SoupX_norm_dataset.Rda')

metadata_short <- metadata[metadata$CellType.Corrected %in% c( "Basal","Cycling Basal", "Suprabasal", "Suprabasal N",
                                                              "Secretory", "Secretory N", "Deuterosome",
                                                              "Multiciliated", "Multiciliated N"),]

markers <- c("MKI67","TOP2A", "CDC20","PCNA","KRT19", "TACSTD2", "KRT5", 
             "BCAM", "LAMB3", "KRT15", "MIR205HG",
             "SLC25A6","MT1X", "CSTA", "SERPINB3","SLC25A5" , "MGST1", "TXN",
             "SCGB1A1", "SLPI")

flt_countTable <- countTable[markers, rownames(metadata_short)]
df_meanExpr <- as.data.frame(matrix(0, ncol = 1+length(unique(metadata_short$CellType.Corrected)), 
                                    nrow = length(markers)))
df_perc <- as.data.frame(matrix(0, ncol = 1 + length(unique(metadata_short$CellType.Corrected)), 
                                nrow = length(markers)))
colnames(df_meanExpr) <- colnames(df_perc) <- c("Gene", unique(metadata_short$CellType.Corrected))
df_meanExpr$Gene <- df_perc$Gene <- rownames(flt_countTable)
for (j in unique(metadata_short$CellType.Corrected)){
  cells_id = rownames(metadata_short[metadata_short$CellType.Corrected == j,])
  df_meanExpr[, j] <- rowMeans(flt_countTable[, cells_id])
  df_perc[, j] <- (rowSums(flt_countTable[, cells_id] > 0) / length(cells_id)) *100
}

meanExpr_long <- gather(df_meanExpr, CellType, Expression, "Suprabasal N":"Multiciliated N")
perc_long <- gather(df_perc, CellType, Percentage, "Suprabasal N":"Multiciliated N")

dotplot_df <- cbind(meanExpr_long, Percentage = perc_long$Percentage)
dotplot_df$CellType <- factor(dotplot_df$CellType, 
                              levels = c("Basal", "Cycling Basal", "Suprabasal", "Suprabasal N",
                                         "Secretory", "Secretory N", "Deuterosome",
                                         "Multiciliated", "Multiciliated N"))

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

# pdf("markers_cycling_cells.pdf", width=7, height=3, useDingbats=FALSE)
# dev.off()

