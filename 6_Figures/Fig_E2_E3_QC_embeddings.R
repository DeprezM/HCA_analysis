## Fig Supp QC

library(ggplot2)
library(tidyr)


setwd("/HCA_analysis/6_Figures/Fig_supp_QC")

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
                     'Endothelial', 'Fibroblast', 'Peicyte', 'Smooth muscle',
                     'B cells', 'Plasma cells', 'LT/NK', 'Mast cells', 'Dendritic', 'Monocyte', 'Macrophage')

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

##########################################################################################
#                                        Plots                                           #
##########################################################################################

### UMAP ---------------------------------------------------------------------------------

ggplot(metadata, aes(x = x, y = y, color = Position))+
  geom_point(size = 0.2) + 
  scale_color_manual(values = position_color) +
  theme_classic() +
  theme(legend.position = "none")

# tiff("umap_position.tiff", width = 4, height = 4, units = "in", res = 800, compression = "jpeg")
# dev.off()

ggplot(metadata, aes(x = x, y = y, color = Method))+
  geom_point(size = 0.2) + 
  scale_color_manual(values = method_color) +
  theme_classic() +
  theme(legend.position = "none")

# tiff("umap_method.tiff", width = 4, height = 4, units = "in", res = 800, compression = "jpeg")
# dev.off()


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Uncorrected UMAP
uncorrected_coord <- read.table(file = "/Data/Uncorrected_embedding.tsv",
                                sep = "\t", header = T, row.names = 1)
colnames(uncorrected_coord)
metadata$unc_umap1 <- uncorrected_coord$umap_1
metadata$unc_umap2 <- uncorrected_coord$umap_2

ggplot(metadata, aes(x = unc_umap1, y = unc_umap2, color = CellType.Corrected_rare))+
  geom_point(size= 0.2) +
  scale_color_manual(values = cell_type_color) +
  theme_classic() +
  theme(legend.position = "none")

# pdf("umap_uncorrected_cell_type.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()

ggplot(metadata, aes(x = unc_umap1, y = unc_umap2, color = Method))+
  geom_point(size= 0.2) +
  scale_color_manual(values = method_color) +
  theme_classic() +
  theme(legend.position = "none")

# pdf("umap_uncorrected_method.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()


ggplot(metadata, aes(x = unc_umap1, y = unc_umap2, color = Position))+
  geom_point(size= 0.2) +
  scale_color_manual(values = position_color) +
  theme_classic() +
  theme(legend.position = "none")

# pdf("umap_uncorrected_position.pdf", width=3, height=3, useDingbats=FALSE)
# dev.off()


### Violin Quality metrics ------------------------------------------------------------------


ggplot(metadata[order(metadata$Position, metadata$Method), ], 
       aes(x = Sample, y = log10(UMI.Count), fill = Position)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = position_color) +
  scale_x_discrete(limits = unique(metadata[order(metadata$Position, metadata$Method), "Sample"])) +
  theme_classic() +
  theme(
    legend.position = "none"
  )
# tiff("UMI_per_sample.tiff", width = 8, height = 4, units = "in", res = 800, compression = "jpeg")
# dev.off()



ggplot(metadata[order(metadata$Position, metadata$Method), ], 
       aes(x = Sample, y = Expressed.Genes, fill = Position)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = position_color) +
  scale_x_discrete(limits = unique(metadata[order(metadata$Position, metadata$Method), "Sample"])) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# tiff("Genes_per_sample.tiff", width = 8, height = 4, units = "in", res = 800, compression = "jpeg")
# dev.off()



### Bar plot - Cell size ---------------------------------------------------------------------

cell_size = data.frame(manip = unique(metadata$Sample),
                       position = substr(unique(metadata$Sample), 11, 13),
                       donor = substr(unique(metadata$Sample), 1, 4),
                       size = c(14.9, 9.46, 7.19,
                                12.6, 8.02, 7.02,
                                12.57,
                                12.96, 9.97, 10.47, 7.12,
                                13.47, 6.25, 9.29, 10.24,
                                13.51, 10.36, 7.92, 6.12,
                                8.52, 8.63, 7.79,
                                9.55, NA, 9.55, 9.55,
                                12.24, 13.11, 12.38, 8.85,
                                11.29, 11.47, 12.69, 12.09, 9.75))

cell_size$position <- factor(cell_size$position, levels = c("Nas", "Pro","Int", "Dis"))

position_color <- c('Nas' = '#ffd966',
                    'Pro' = '#e69138',
                    'Int' = '#e06666',
                    'Dis' = '#85200c')

# mean(cell_size[cell_size$position == "Dis",]$size, na.rm = T)
# sd(cell_size[cell_size$position == "Dis",]$size, na.rm = T)

ggplot(cell_size, aes(x = position, y = size, fill = position))+
  geom_boxplot()+
  geom_jitter(size = 2)+
  scale_fill_manual(values = position_color) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_blank())

pdf("sample_cell_size.pdf", width=4, height=3, useDingbats=FALSE)
dev.off()



### Statistical tests QC --------------------------------------------------------------------------

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# Nb of genes
t.test(metadata[metadata$Method == "Biopsy" & metadata$Position == "Nasal",]$Expressed.Genes,
       metadata[metadata$Method == "Brushing" & metadata$Position == "Nasal",]$Expressed.Genes)

boxplot(metadata[metadata$Method == "Biopsy" & metadata$Position == "Nasal",]$Expressed.Genes,
        metadata[metadata$Method == "Brushing" & metadata$Position == "Nasal",]$Expressed.Genes)


t.test(metadata[metadata$Method == "Biopsy" & metadata$Position == "Nasal",]$Expressed.Genes,
       metadata[metadata$Method == "Biopsy" & metadata$Position == "Proximal",]$Expressed.Genes)

t.test(metadata[metadata$Method == "Biopsy" & metadata$Position == "Proximal",]$Expressed.Genes,
       metadata[metadata$Method == "Biopsy" & metadata$Position == "Intermediate",]$Expressed.Genes)

t.test(metadata[metadata$Method == "Biopsy" & metadata$Position == "Proximal" & 
                  metadata$Donor %in% c("D326"),]$Expressed.Genes,
       metadata[metadata$Method == "Biopsy" & metadata$Position == "Intermediate" & 
                  metadata$Donor %in% c("D326"),]$Expressed.Genes)

boxplot(metadata[metadata$Method == "Biopsy" & metadata$Position == "Proximal",]$Expressed.Genes,
        metadata[metadata$Method == "Biopsy" & metadata$Position == "Intermediate",]$Expressed.Genes)



t.test(metadata[metadata$Method == "Brushing" & metadata$Position == "Distal",]$Expressed.Genes,
       metadata[metadata$Method == "Brushing" & metadata$Position == "Nasal",]$Expressed.Genes)

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# Nb of UMI

t.test(metadata[metadata$Method == "Biopsy" & metadata$Position == "Nasal",]$UMI.Count,
       metadata[metadata$Method == "Brushing" & metadata$Position == "Nasal",]$UMI.Count)

boxplot(metadata[metadata$Method == "Biopsy" & metadata$Position == "Nasal",]$UMI.Count,
        metadata[metadata$Method == "Brushing" & metadata$Position == "Nasal",]$UMI.Count)


t.test(metadata[metadata$Method == "Biopsy" & metadata$Position == "Nasal",]$UMI.Count,
       metadata[metadata$Method == "Biopsy" & metadata$Position == "Proximal",]$UMI.Count)

t.test(metadata[metadata$Method == "Biopsy" & metadata$Position == "Proximal",]$UMI.Count,
       metadata[metadata$Method == "Biopsy" & metadata$Position == "Intermediate",]$UMI.Count)

t.test(metadata[metadata$Method == "Biopsy" & metadata$Position == "Proximal" & 
                  metadata$Donor %in% c("D326"),]$UMI.Count,
       metadata[metadata$Method == "Biopsy" & metadata$Position == "Intermediate" & 
                  metadata$Donor %in% c("D326"),]$UMI.Count)

boxplot(metadata[metadata$Method == "Biopsy" & metadata$Position == "Proximal",]$UMI.Count,
        metadata[metadata$Method == "Biopsy" & metadata$Position == "Intermediate",]$UMI.Count)



t.test(metadata[metadata$Method == "Brushing" & metadata$Position == "Distal",]$UMI.Count,
       metadata[metadata$Method == "Brushing" & metadata$Position == "Nasal",]$UMI.Count)



# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #
# Cell size

wilcox.test(cell_size[cell_size$position == "Nas",]$size,
            cell_size[cell_size$position == "Pro",]$size )

wilcox.test(cell_size[cell_size$position == "Int",]$size,
            cell_size[cell_size$position == "Dis",]$size )

wilcox.test(cell_size[cell_size$position == "Nas",]$size,
            cell_size[cell_size$position == "Dis",]$size )









