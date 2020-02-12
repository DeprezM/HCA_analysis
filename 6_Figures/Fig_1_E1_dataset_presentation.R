library(ggplot2)
library(tidyr)

setwd("/HCA_analysis/6_Figures/Figure_complete_dataset")

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


summary_table <- as.data.frame(table(metadata[order(metadata$Position, metadata$Method),]$Sample, metadata[order(metadata$Position, metadata$Method),]$CellType.Corrected))
head(summary_table)
summary_table <- spread(summary_table, Var2, Freq)
rownames(summary_table) <- summary_table$Var1

summary_table <- summary_table[unique(metadata[order(metadata$Position, metadata$Method),]$Sample),]
# write.table(summary_table, file = "Nb_cell_sample_celltype.tsv", sep = "\t", quote = F, row.names = F)


##########################################################################################
#                                        Plots                                           #
##########################################################################################

### UMAP ---------------------------------------------------------------------------------

ggplot(metadata, aes(x = x, y = y, color = CellType.Corrected_cycling))+
  geom_point(size = 0.2) + 
  scale_color_manual(values = cell_type_color) +
  theme_classic() +
  theme(legend.position = "none")

# pdf("umap_cell_type.pdf", width=4, height=4, useDingbats=FALSE)
# dev.off()


ggplot(metadata, aes(x = x, y = y, color = Donor))+
  geom_point(size = 0.2, alpha = 0.5) + 
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~Donor)

# pdf("umap_donor.pdf", width=4, height=4, useDingbats=FALSE)
# pdf("umap_donor_large.pdf", width=15, height=15, useDingbats=FALSE)
# dev.off()

ggplot(metadata, aes(x = x, y = y, color = Sample))+
  geom_point(size = 0.2, alpha = 0.5) + 
  theme_classic() +
  theme(legend.position = "none")+
  facet_wrap(~Sample)

# pdf("umap_sample.pdf", width=4, height=4, useDingbats=FALSE)
# pdf("umap_sample_large.pdf", width=15, height=15, useDingbats=FALSE)
# dev.off()

### Proportion ---------------------------------------------------------------------------

sampletype_order <- c('Nasal\nBrushing', 'Nasal\nBiopsy', 
                      'Proximal\nBiopsy', 'Intermediate\nBiopsy',
                      'Distal\nBrushing')

metadata$CellType.Corrected_cycling <- factor(metadata$CellType.Corrected_cycling, levels = cell_type_order)
metadata$Position <- factor(metadata$Position, levels = names(position_color))
metadata$SampleType = paste0(metadata$Position, '\n', metadata$Method)
metadata$SampleType = factor(metadata$SampleType, levels = sampletype_order)


df = as.data.frame.matrix(table(metadata$SampleType, metadata$CellType.Corrected_cycling))
df = as.data.frame(df / rowSums(df))
df$Sample <- factor(rownames(df), 
                    levels = unique(metadata[order(metadata$SampleType, metadata$Position, metadata$Method), "Sample"]))
df$Sample <- factor(rownames(df), 
                    levels = sampletype_order)
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

# pdf("proportion_cell_type.pdf", width = 7, height = 7, useDingbats=FALSE)
# dev.off()


df = as.data.frame(table(metadata$CellType.Corrected_cycling))
df$Percentage = df$Freq / sum(df$Freq)
df$Sample <- "Dataset"

ggplot(df, aes(x = "", y = Percentage, fill = Var1)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0, direction = -1) +
  scale_fill_manual(values = cell_type_color) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none" 
  )

# pdf("proportion_cell_type_in_dataset.pdf", width = 4, height = 4, useDingbats=FALSE)
# dev.off()


df = as.data.frame.matrix(table(metadata$Sample, metadata$CellType.Corrected_rare))
df = as.data.frame(df / rowSums(df))
df$Sample <- factor(rownames(df), 
                    levels = unique(metadata[order(metadata$SampleType, metadata$Position, metadata$Method), "Sample"]))
df_long = gather(df, cell_type, percentage, 'B cells':'Suprabasal N')
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

# pdf("proportion_cell_type_per_sample.pdf", width = 8, height = 4, useDingbats=FALSE)
# dev.off()


### Heatmap ---------------------------------------------------------------------------------

load('/Data/soupX_norm_dataset.Rda')

markers <- c("KRT5", "TP63", "DLK2", # Basal
             "SERPINB4", "KRT19", "MGST1", # Suprabasal "SERPINB4", "NTS", "KRT19",
             "FABP5", "TFCP2L1", "S100A9",  # Suprabasal N "SAT1", "SGK1",
             "SCGB1A1", "TSPAN8", "MUC5AC", # Secretory
             "MUC1", "S100P", "PSCA", # Secretory N
             "CDC20Bshort", "HES6", "CCDC67", # Deuterosomal
             "FOXJ1", "PIFO", "RYR3", # Multiciliated
             "TPRXL", "PALLD", "CYP24A1", # Multiciliated N
             "LTF", "LYZ", "PIP", # Serous
             "AZGP1", "MUC5B", "BPIFB2",  # SMG Goblet
             "AQP1", "GNG11", "DARC", # Endothelial
             "FBLN1", "DCN", "C1R", # Fibroblast
             "DES", "CNN1", "ACTA2",  # Smooth muscle
             "MYL9", "RERGL", "PDGFRB", # Pericyte
             "TYROBP", "APOC1", "C1QA", # Macrophage
             "SDS", "C15orf48", "FNIP2", # Monocyte
             "CST3", "MS4A6A", "HLA-DPB1", # Dendritic
             "TPSAB1", "CPA3", "HPGDS", # Mast cells
             "IL32", "CD3D", "CCL5", # LT/NK
             "LTB", "MS4A1", "CD79A", # B cells
             "IGJ", "SSR4", "MZB1" # Plasma cells
)

markers_flt_count <- countTable[markers, ]
colnames(markers_flt_count) <- gsub(".", "-", colnames(countTable), fixed = T)
colnames(countTable) <- gsub(".", "-", colnames(countTable), fixed = T)

markers_count <- matrix(0, nrow = length(markers), ncol = length(cell_type_order) - 8)
colnames(markers_count) <- c("Basal", "Suprabasal", "Suprabasal N", "Secretory", 
                             "Secretory N", "Deuterosomal", "Multiciliated", "Multiciliated N",
                             "Serous", "SMG Goblet", "Endothelial", "Fibroblast",
                             "Smooth muscle", "Pericyte","Macrophage", "Monocyte", "Dendritic",  
                             "Mast cells", "LT/NK", "B cells", "Plasma cells")

rownames(markers_count) <- markers

for (i in 1:ncol(markers_count)){
  print(colnames(markers_count)[i])
  cell_names <- rownames(metadata[metadata$CellType.Corrected_cycling == colnames(markers_count)[i], ])
  markers_count[, colnames(markers_count)[i]] <- rowMeans(markers_flt_count[, cell_names])
}


markers_count <- flexible_normalization(markers_count, by_row = T)
markers_count[markers_count > 2] <- 2
markers_count[markers_count < -1.5] <- -1.5

annotation <- data.frame(CellType = factor(colnames(markers_count), levels = cell_type_order))
rownames(annotation) <- colnames(markers_count)
ann_colors <- list(CellType = cell_type_color)

pp <- pheatmap::pheatmap(markers_count, cluster_rows = F, cluster_cols = F, 
                         annotation_col = annotation, annotation_colors = ann_colors,     
                         angle_col = 90, annotation_legend = F)

# save_pheatmap_pdf(pp, "Heatmap_markers_1.pdf", width = 7, height = 8)



### Violin gene expression ---------------------------------------------------------------------


metadata$gene <- as.numeric(countTable["BPIFA2",])
ggplot(metadata[order(metadata$Position, metadata$Method), ], 
       aes(x = Sample, y = gene)) +
  geom_violin(scale = "width", fill = "#84C553") + 
  scale_x_discrete(limits = unique(metadata[order(metadata$Position, metadata$Method), "Sample"])) +
  geom_jitter(size = 0.2) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))

# pdf("violin_gene_sample.pdf", width=8, height=4, useDingbats=FALSE)
# dev.off()

