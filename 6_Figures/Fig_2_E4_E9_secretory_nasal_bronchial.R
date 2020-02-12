# Differential analysis

library(ggplot2)
library(tidyr)
library(gridExtra)
library(corrplot)
library(VennDiagram)
library(fgsea)
library(gage)
library(biomaRt)
library(ggrepel)

suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggrepel))

setwd("/HCA_analysis/6_Figures/Figure_nasal_bronchial")

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


panel.cor_simple <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- cor(x, y) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * abs(r)) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
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



##########################################################################################
#                                        Plots                                           #
##########################################################################################

### MA-plot ------------------------------------------------------------------------------

topTable <- read.table("DA_bulk_secretory.csv", 
                       sep = "\t", header = T)
# topTable <- read.table("DA_bulk_suprabasal.csv", 
#                        sep = "\t", header = T)
# topTable <- read.table("DA_bulk_multiciliated.csv", 
#                        sep = "\t", header = T)
topTable$toplot <- rep(0,nrow(topTable))
topTable$toplot[topTable$FDR<0.05 & topTable$logFC > 1.5 & topTable$logCPM > 1.5 ] <- 1
topTable$toplot[topTable$FDR<0.05 & topTable$logFC < -1.5 & topTable$logCPM > 1.5 ] <- -1


topTable$towrite[topTable$FDR<0.05 & abs(topTable$logFC) > 3 & abs(topTable$logCPM) > 4 ] <- 1
table(topTable$towrite)
rownames(topTable) <- topTable$Gene

rownames(topTable[topTable$towrite == 1,])

topTable$towrite <- rep(0,nrow(topTable))
topTable[c("ITLN1","SCGB1A1","CTD-2531D15.4","SPRR1A","DUOX2","MUC5B","TMEM45A","RP11-89K21.1","RP11-355F16.1",
           "SIX3","SIX3-AS1","DUOXA2","PAX7","CLCA2","APOC1","SCGB3A1","TFF3","MS4A8",
           "C16orf89","CYP1B1","CLEC7A","SERPINB7","SYT8","CD36","ATP10B","LYNX1","SCNN1G" ,
           "SERPINF1","BEND5","BPIFB1","AGR3","KLK11","CXCL6","SYT7","VIM","CLDN10" ,
           "TMEM241","C1orf168","ADRA2A","GSTA2","GDF15","CHP2","XDH","GALNT6",
           "GPX2","FRMD4A", "CA12","S100A4","MT1G","S100A8","CCDC64","ASRGL1","BPIFA2","PI3","MUC4"
), "towrite"] <- 1

table(topTable$toplot)

gg <- ggplot(data=topTable, aes(x=logCPM,
                                y=logFC)) +  
  geom_point(color="grey60", fill = "grey60", size=1, shape = 21) +
  geom_point(data=subset(topTable,toplot==1),shape = 21, fill = "red", colour = "red4", size=1.3, alpha = 0.8) +
  geom_point(data=subset(topTable,toplot==-1),shape = 21, fill = "blue", colour = "blue4", size=1.3, alpha = 0.8) +
  theme_classic() +
  geom_point(data=subset(topTable,towrite==1), shape = 21, colour = "black", size=2) +
  geom_text_repel(data=subset(topTable,towrite==1) ,aes(label=Gene), size=3, 
                  show.legend=FALSE, force=2)
gg 

pdf("MA_plot_secretory.pdf", width = 5, height = 5, useDingbats=FALSE)
# pdf("MA_plot_suprabasal.pdf", width = 5, height = 5, useDingbats=FALSE)
# pdf("MA_plot_multiciliated.pdf", width = 5, height = 5, useDingbats=FALSE)
dev.off()


### Violin plot ------------------------------------------------------------------------------

nasal_genes <- c("MUC4", "PI3", "DUOX2", "PAX7", "SIX3", "FOXG1", "LYNX1",
                 "S100A4", "CEACAM5", "LYPD2", "BPIFA2", "STATH", "LY6D",
                 "SDCBP2", "ASRGL1")

bronchial_genes <- c("CLCA2",  "MUC5B", "SCGB1A1",  "SCGB3A1", "TMEM45A", "C16orf89", "TFF3",   "AGR3","KLK11",
                     "SERPINB7", "SERPINF1", "CLDN10",  "CXCL6","FOLR1",
                     "GSTA2")


gene <- nasal_genes
# gene <- bronchial_genes

metadata_genes <- cbind(metadata, t(countTable[gene,]))
meta_subset <- metadata_genes[metadata_genes$CellType.Corrected %in% c("Secretory", "Secretory N"),]
meta_subset <- meta_subset[, c("SampleType", gene)]
df_long <- gather(meta_subset, Gene, Value, gene)
df_long$Gene <- factor(df_long$Gene, levels = gene)
head(df_long)

ggplot(df_long, 
       aes(x = SampleType, y = Value)) +
  geom_violin(scale = "width", fill = "#84C553") + 
  scale_x_discrete(limits = unique(sampletype_order)) +
  scale_y_continuous(position = "right", breaks = c(min(df_long$Value), max(df_long$Value)))+
  facet_wrap(~Gene, ncol = 1, strip.position="left", scales = "free_y") +
  theme_classic()+
  theme(
    strip.background = element_rect(colour="white", fill="white"),
    axis.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 180),
    # axis.text.y = element_text(size = 15),
    axis.title = element_blank())

# pdf(paste0("Violin_sampletype_secretory_wrap_nasal.pdf"), width = 4, height = 9, useDingbats=FALSE)
# dev.off()


### Venn Diagram ------------------------------------------------------------------------------


tt_secretory <- read.table("DA_bulk_secretory.tsv", 
                           sep = "\t", header = T)
tt_suprabasal <- read.table("DA_bulk_suprabasal.tsv", 
                            sep = "\t", header = T)
tt_multiciliated <- read.table("DA_bulk_multiciliated.tsv", 
                               sep = "\t", header = T)

## Nasal genes
SCs_are <- as.character(tt_secretory[tt_secretory$FDR<0.05 & tt_secretory$logFC > 3 & tt_secretory$logCPM > 2 , "Gene"])
MCCs_are <- as.character(tt_multiciliated[tt_multiciliated$FDR<0.05 & tt_multiciliated$logFC > 3 & tt_multiciliated$logCPM > 2 , "Gene"])
SBs_are <-  as.character(tt_suprabasal[tt_suprabasal$FDR<0.05 & tt_suprabasal$logFC > 3 & tt_suprabasal$logCPM > 2 , "Gene"])

length(unique(c(SCs_are, SBs_are, MCCs_are)))
draw.triple.venn(area1 = length(SCs_are), area2 = length(MCCs_are), area3 = length(SBs_are),
                 n12 = length(intersect(SCs_are, MCCs_are)), 
                 n23 = length(intersect(MCCs_are, SBs_are)), 
                 n13 = length(intersect(SCs_are, SBs_are)), 
                 n123 = length(intersect(SCs_are, intersect(MCCs_are, SBs_are))), 
                 category = c("Secretory N", "Multiciliated N", "Suprabasal N"), lty = "blank", 
                 fill = c("#a9c653", "#7e98ce", "#e0c96c"))

# pdf("Venn_up_Nasal.pdf", width = 5, height = 5, useDingbats=FALSE)
# dev.off()

## Bronchial genes
SCs_are <- as.character(tt_secretory[tt_secretory$FDR<0.05 & tt_secretory$logFC < -3 & tt_secretory$logCPM > 2 , "Gene"])
MCCs_are <- as.character(tt_multiciliated[tt_multiciliated$FDR<0.05 & tt_multiciliated$logFC < -3 & tt_multiciliated$logCPM > 2 , "Gene"])
SBs_are <-  as.character(tt_suprabasal[tt_suprabasal$FDR<0.05 & tt_suprabasal$logFC < -3 & tt_suprabasal$logCPM > 2 , "Gene"])


draw.triple.venn(area1 = length(SCs_are), area2 = length(MCCs_are), area3 = length(SBs_are),
                 n12 = length(intersect(SCs_are, MCCs_are)), 
                 n23 = length(intersect(MCCs_are, SBs_are)), 
                 n13 = length(intersect(SCs_are, SBs_are)), 
                 n123 = length(intersect(SCs_are, intersect(MCCs_are, SBs_are))), 
                 category = c("Secretory", "Multiciliated", "Suprabasal"), lty = "blank", 
                 fill = c("#53c653", "#466cb9", "#FCCC0A"))

# pdf("Venn_up_Bronchial.pdf", width = 5, height = 5, useDingbats=FALSE)
# dev.off()



### Bar Plot GSEA ------------------------------------------------------------------------------

df_GO_bp <- read.table("secretory_NASAL_fgseaRes_GO_bp.csv", 
                      sep = "\t", header = T)

# secretory_BRONCHIAL_fgseaRes_GO_bp.csv
# suprabasal_NASAL_fgseaRes_GO_bp.csv
# suprabasal_BRONCHIAL_fgseaRes_GO_bp.csv
# multiciliated_NASAL_fgseaRes_GO_bp.csv
# multiciliated_BRONCHIAL_fgseaRes_GO_bp.csv
# 

# Secretory Bronchial
paths <- c("GO_DEFENSE_RESPONSE",
           "GO_RESPONSE_TO_INORGANIC_SUBSTANCE",
           "GO_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
           "GO_INNATE_IMMUNE_RESPONSE",
           "GO_REGULATION_OF_RESPONSE_TO_WOUNDING")

pathToPlot <- df_GO_bp[df_GO_bp$pathway %in% paths,]
pathToPlot$pathway <- factor(pathToPlot$pathway, levels = paths)

ggplot(pathToPlot, aes(y = -log(pval), x = pathway)) +
  geom_bar(stat = "identity", fill = "lightcyan4") +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        # axis.text = element_text(size = 13, face = "bold"), 
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

pdf("Gsea_secretory_bronchial_pval.pdf", width =3.5, height = 2, useDingbats=FALSE)
dev.off()

ggplot(pathToPlot, 
       aes(y = -Fraction, x = pathway)) +
  geom_bar(stat = "identity", fill = "dodgerblue4") +
  ylim(-100, 0) +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text = element_text(size = 13, face = "bold"),
        legend.position = "none")
pdf("Gsea_secretory_bronchial_fraction.pdf", width =5, height = 2, useDingbats=FALSE)
dev.off()



# Secretory NASAL
paths <- c("GO_EPIDERMAL_CELL_DIFFERENTIATION",
           "GO_REGIONALIZATION",
           "GO_SENSORY_PERCEPTION",
           "GO_EPITHELIAL_CELL_DIFFERENTIATION",
           "GO_CELL_MOTILITY")

pathToPlot <- df_GO_bp[df_GO_bp$pathway %in% paths,]
pathToPlot$pathway <- factor(pathToPlot$pathway, levels = paths)

ggplot(pathToPlot, aes(y = -log(pval), x = pathway)) +
  geom_bar(stat = "identity", fill = "lightcyan4") +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        # axis.text = element_text(size = 13, face = "bold"), 
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

pdf("Gsea_secretory_nasal_pval.pdf", width =3.5, height = 2, useDingbats=FALSE)
dev.off()

ggplot(pathToPlot, 
       aes(y = -Fraction, x = pathway)) +
  geom_bar(stat = "identity", fill = "dodgerblue4") +
  ylim(-100, 0) +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        # axis.text = element_text(size = 13, face = "bold"),
        legend.position = "none")
pdf("Gsea_secretory_nasal_fraction.pdf", width =5, height = 2, useDingbats=FALSE)
dev.off()




### Heatmap correlation cell types -----------------------------------------------------------

metadata$nasal_celltype <- as.character(metadata$CellType.Corrected)

celltype <- c("Basal", "Suprabasal", "Secretory", "Multiciliated", "Deuterosomal", 
              "Suprabasal N", "Secretory N", "Multiciliated N", "Serous", "SMG Goblet")
meta <- metadata[metadata$nasal_celltype %in% celltype, ]

celltype_B <- c("Basal", "Suprabasal", "Secretory", "Multiciliated", "Deuterosomal","Serous", "SMG Goblet")
meta$position_cell_type <- meta$nasal_celltype
meta$position_cell_type[meta$nasal_celltype %in% celltype_B & 
                          meta$Position == "Nasal"] <- 
  paste0(meta$position_cell_type[meta$nasal_celltype %in% celltype_B & meta$Position == "Nasal"], " N")
meta$position_cell_type[meta$nasal_celltype %in% celltype_B & 
                          meta$Position == "Proximal"] <- 
  paste0(meta$position_cell_type[meta$nasal_celltype %in% celltype_B & meta$Position == "Proximal"], " P")
meta$position_cell_type[meta$nasal_celltype %in% celltype_B & meta$Position == "Intermediate"] <- 
  paste0(meta$position_cell_type[ meta$nasal_celltype %in% celltype_B & meta$Position == "Intermediate"], " I")
meta$position_cell_type[meta$nasal_celltype %in% celltype_B & meta$Position == "Distal"] <- 
  paste0(meta$position_cell_type[meta$nasal_celltype %in% celltype_B & meta$Position == "Distal"], " D")

cells_id = rownames(meta)

for (i in cells_id){

  meta[i, ]$position_cell_type = paste0(unlist(strsplit( meta[i, ]$position_cell_type, " "))[1], " ",
                                        substr(meta[i, ]$Position, 1, 1))
}

table(meta$position_cell_type)
meta$manip_cell_type <- paste(meta$position_cell_type, meta$Sample)
celltype = unique(meta$manip_cell_type)
rownames(meta) <- gsub("-",".", rownames(meta))

load('/Data/scranNorm_dataset.Rda')

norm = countTable[, rownames(meta)]
dim(norm)

cmed = c()
kept_sample <- c()
for (i in 1:length(celltype)) {
  print(celltype[i])
  if (sum(meta$manip_cell_type == celltype[i]) > 20) {
    cmed = cbind(cmed, rowMeans(norm[, meta$manip_cell_type==celltype[i]]))
    kept_sample <- c(kept_sample, celltype[i])
  }
}
colnames(cmed)=kept_sample
ct = cor(cmed)

unlist_el <- unlist(strsplit(colnames(ct), " "))
annotation <- data.frame(cell_type = unlist_el[seq(1,length(unlist_el), 3)],
                         position = unlist_el[seq(2,length(unlist_el), 3)],
                         donor = unlist_el[seq(3,length(unlist_el), 3)],
                         stringsAsFactors = F)

rownames(annotation) <- colnames(ct)
annotation["Secretory N", "cell_type"] <- "Secretory N"
annotation["Suprabasal N", "cell_type"] <- "Suprabasal N"


position_color <- c('N' = '#ffd966',
                    'P' = '#e69138',
                    'I' = '#e06666',
                    'D' = '#85200c')

cell_type_color["SMG"] <- cell_type_color["SMG Goblet"]
anno_col <- list(cell_type = cell_type_color, position = position_color)

pheatmap(ct, labels_row = F)
pp <- pheatmap(ct, main="",labels_row = F, annotation_legend = F,
               annotation_col = annotation, annotation_colors = anno_col)

# save_pheatmap_pdf(pp, "Heatmap_correlation_epithelial_cell_type_per_sample.pdf", width = 8.5, height = 9)


### Heatmap Score TF -----------------------------------------------------------


load("/Data/TFscore_matrix_.Rda")
# load("/Data/Secretory_TFTopTable.Rda")

tf_matrix <- as.matrix(TF_matrix@assays$RNA@counts)

markers <- c(rownames(A[order(A$avg_logFC, decreasing = T), ])[1:20],"AHR",
             rownames(A[order(A$avg_logFC, decreasing = F), ])[1:11])

markers_flt_count <- tf_matrix[markers, ]
colnames(markers_flt_count) <- gsub(".", "-", colnames(countTable), fixed = T)

markers_count <- matrix(0, nrow = nrow(markers_flt_count), ncol = 10)
colnames(markers_count) <- c("Basal", "Suprabasal", "Suprabasal N", "Secretory", 
                             "Secretory N", "Deuterosomal", "Multiciliated", "Multiciliated N",
                             "Serous", "SMG Goblet")

rownames(markers_count) <- markers

for (i in 1:ncol(markers_count)){
  print(colnames(markers_count)[i])
  cell_names <- rownames(metadata[metadata$CellType.Corrected == colnames(markers_count)[i], ])
  markers_count[, colnames(markers_count)[i]] <- rowMedians(markers_flt_count[, cell_names])
}

markers_count[markers_count > 1] <- 1

annotation <- data.frame(CellType = factor(colnames(markers_count), levels = cell_type_order))
rownames(annotation) <- colnames(markers_count)
ann_colors <- list(CellType = cell_type_color)

pp <- pheatmap::pheatmap(markers_count, cluster_rows = F, cluster_cols = F, 
                         annotation_col = annotation, annotation_colors = ann_colors,     
                         angle_col = 90, annotation_legend = F)

# save_pheatmap_pdf(pp, "Diff_secretory_TF.pdf", width = 7, height = 8)

### UMAP focus secretory  -----------------------------------------------------------

load('/Data/SoupX_norm_dataset.Rda')
secretory_umap <- read.csv(file = "/Data/Secretory_UMAP.tsv",
                           sep = "\t", stringsAsFactors = F, header = T, row.names = 1)

secretoryN_umap <- read.csv(file = "/Data/SecretoryN_UMAP.tsv",
                            sep = "\t", stringsAsFactors = F, header = T, row.names = 1)


gene <- c("SCGB1A1", "SCGB3A1", "MUC5AC", "MUC5B")
secretory_umap <- cbind(secretory_umap, t(countTable[gene,rownames(secretory_umap)]))
long <- gather(secretory_umap, Gene, Value, "SCGB1A1":"MUC5B")

ggplot(long, aes(x = X0, y = X1, color = Value))+ 
  geom_point(size = 0.3)+
  scale_color_gradient(low = "grey70", high = "red4") +
  facet_wrap(~Gene, scales = "free") + 
  theme_classic() + theme(legend.position = "none")

# pdf("Umap_secretory_gradient.pdf", width = 4, height = 4, useDingbats=FALSE)
# dev.off()

secretoryN_umap <- cbind(secretoryN_umap, t(countTable[gene,rownames(secretoryN_umap)]))
long <- gather(secretoryN_umap, Gene, Value, "SCGB1A1":"MUC5B")

ggplot(long, aes(x = X0, y = X1, color = Value))+ 
  geom_point(size = 0.3)+
  scale_color_gradient(low = "grey70", high = "red4") +
  facet_wrap(~Gene, scales = "free") + 
  theme_classic() + theme(legend.position = "none")


# pdf("Umap_secretoryN_gradient.pdf", width = 4, height = 4, useDingbats=FALSE)
# dev.off()


