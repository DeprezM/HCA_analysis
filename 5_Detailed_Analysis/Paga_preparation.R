## PAGA Analysis

library(ggplot2)
library(tidyr)
library(scales)
library(forcats)
library(matrixStats)
library(fgsea)

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



## Load Rare cells metadata -------------------------------------------------------------

meta_rare <- metadata[metadata$CellType.Corrected_rare == "Rare cells",]

rare_coords <- read.table(file = "/Data/Focus_rare_cells_umap_coord.tsv", 
                          sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
meta_rare$x_rareCoord <- rare_coords[rownames(meta_rare), 'X0']
meta_rare$y_rareCoord <- rare_coords[rownames(meta_rare), 'X1']
# AT1 and AT2 cells are low abundant cells in this dataset but they are not so-called 'rare cells'
# meta_rare = meta_rare[!meta_rare$CellType.Corrected %in% c("AT1", "AT2"),]



##########################################################################################
#                                 Subset cells for PAGA                                  #
##########################################################################################


cell_type_sub <- c("Cycling Basal", "Basal", "Suprabasal", "Secretory", "Deuterosomal", "Multiciliated")
cell_type_subAll <- c("Cycling Basal", "Basal", "Deuterosomal", "Suprabasal N", "Secretory N", "Multiciliated N")

sub <- c()
for (i in 1:length(cell_type_sub)) {
  if (sum(metadata$CellType.Corrected_rare == cell_type_sub[i] &
          metadata$Position != "Nasal" &
          metadata$Sample != "D339_Biop_Int1") > 500){
    nbCells = 500
  } else {
    nbCells = sum(metadata$CellType.Corrected_rare == cell_type_sub[i] & metadata$Position != "Nasal" &
                    metadata$Sample != "D339_Biop_Int1")
  }
  
  sub <- c(sub, sample(rownames(metadata[metadata$CellType.Corrected_rare == cell_type_sub[i] &
                                           metadata$Position != "Nasal" &
                                           metadata$Sample != "D339_Biop_Int1", ]),  nbCells, replace = F) )
}
sub <- c(sub, rownames(metadata[metadata$Position != "Nasal" & metadata$CellType.Corrected %in% "Mucous Multiciliated cells", ]))

subAll <- c()
for (i in 1:length(cell_type_sub)) {
  if (sum(metadata$CellType.Corrected_rare == cell_type_sub[i] &
          metadata$Position == "Nasal" &
          metadata$Sample != "D322_Biop_Nas1") > 500){
    nbCells = 500
  } else {
    nbCells = sum(metadata$CellType.Corrected_rare == cell_type_sub[i] & metadata$Position == "Nasal"  &
                    metadata$Sample != "D322_Biop_Nas1")
  }
  subAll <- c(subAll, sample(rownames(metadata[metadata$CellType.Corrected_rare == cell_type_subAll[i]  &
                                                 metadata$Position == "Nasal" &
                                                 metadata$Sample != "D322_Biop_Nas1", ]),  nbCells, replace = F) )
}
subAll <- c(subAll, rownames(metadata[metadata$Position == "Nasal" & metadata$CellType.Corrected %in% "Mucous Multiciliated cells", ]))


write.table(sub, file = "/data/deprez_data/HCA/Analysis/FullDataset_v6/PAGA_subset_Bronchial.txt",
            quote = F, row.names = F)
write.table(subAll, file = "/data/deprez_data/HCA/Analysis/FullDataset_v6/PAGA_subset_Nasal.txt",
            quote = F, row.names = F)

table(metadata[sub,]$CellType.Corrected)
table(metadata[sub,]$CellType.Corrected_rare)
