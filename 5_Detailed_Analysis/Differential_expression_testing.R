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

load('/Data/raw_exprMatrix.Rda')
# load('/Data/soupX_raw_dataset.Rda')

dim(countTable)

rownames(metadata)[1]
colnames(countTable)[1]
colnames(countTable) <- gsub(".", "-", colnames(countTable), fixed = T)
countTable <-countTable[, rownames(metadata)]

##########################################################################################
#                Differential expression testing : Secretory cells                       #
##########################################################################################

## Filter out Ribo genes and MT genes
RB_genes = read.delim("/data/deprez_data/HCA/PeerLab_analysis/RB_genes", header = F, stringsAsFactors = F)
sbt_count = countTable[!(rownames(countTable) %in% 
                           c(grep("^MT-", rownames(countTable), value = T), RB_genes[,"V1"])),]

# Estimate the number of concerned cells
table(metadata$CellType.Corrected)
table(metadata$CellType.Corrected, metadata$Sample)
metadata$annotation <- metadata$CellType.Corrected

## Subset dataset ------------------------------------------------------------------------

# Remove sample with very different cell type composition and thus very different gene expression background
sbt_count <- sbt_count[, rownames(metadata[metadata$annotation %in% c("Secretory", "Secretory N") &
                                             !metadata$Sample %in% c("D339_Brus_Dis1", "D363_Brus_Dis1"), ])]
sbt_meta <- metadata[rownames(metadata[metadata$annotation %in% c("Secretory", "Secretory N") &
                                         !metadata$Sample %in% c("D339_Brus_Dis1", "D363_Brus_Dis1"), ]), ]
table(sbt_meta$annotation, sbt_meta$Sample)
table(sbt_meta$annotation)

# Dispatch roughly the same number of cells across the future 'mini-bulks'
# To create permanent cell selection (avoid random) set a seed parameter 
sbt_meta$bulk_id <- ""
sbt_meta[sbt_meta$annotation == "Secretory", "bulk_id"] <- paste0("bulk_SC_", sample(
  rep(1:9, 1000), sum(sbt_meta$annotation == "Secretory"), replace=FALSE))
sbt_meta[sbt_meta$annotation == "Secretory N", "bulk_id"] <- paste0("bulk_SCN_", sample(
  rep(1:10, 1000), sum(sbt_meta$annotation == "Secretory N"), replace=FALSE))

# with soupX data
# sbt_meta[sbt_meta$annotation == "Secretory N", "bulk_id"] <- paste0("bulk_SCN_", sample(
#   rep(1:8, 1500), sum(sbt_meta$annotation == "Secretory N"), replace=FALSE))

# Create the bulk samples
bulk_samples <- matrix(data = 0, nrow = nrow(sbt_count),
                       ncol = length(unique(sbt_meta$bulk_id)))
colnames(bulk_samples) <- unique(sbt_meta$bulk_id)

for (sample in unique(sbt_meta$bulk_id)){
  print(sample)
  ix <- rownames(sbt_meta[sbt_meta$bulk_id == sample,])
  print(length(ix))
  bulk_samples[, sample] <- rowSums(sbt_count[, ix])
}
rownames(bulk_samples) <- rownames(sbt_count)


## Prepare summary file - design for statistical testing
serie <- data.frame(Sample=colnames(bulk_samples),
                    type = unlist(strsplit(colnames(bulk_samples), "_"))[seq(2, ncol(bulk_samples)*3, 3)])

## Check for equivalent library size
serie$Counts <- colSums(bulk_samples, na.rm = T)
gg <-ggplot(serie, aes(Sample, Counts)) + 
  geom_bar(stat = "identity", aes(fill=type)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
gg

# pdf("library_size_bulk_secretory.pdf", width=7, height=4, useDingbats=FALSE)
# pdf("library_size_bulk_secretory_soupX.pdf", width=7, height=4, useDingbats=FALSE)
# dev.off()

## Correlation between cluster
mycols=brewer.pal(9,"Spectral")
names(mycols) <- unique(serie$type)
annotation_col = data.frame(Type = serie$type,
                            Position = serie$Sample)
rownames(annotation_col) = serie$Sample
ann_colors = list(Type=mycols)

log_bulk_sample <- log(bulk_samples + 1)
ct <- cor(log_bulk_sample)
pp<- pheatmap(ct, show_rownames = F, show_colnames = T, 
         annotation_col = annotation_col, annotation_colors = ann_colors)
# save_pheatmap_pdf(pp, "Correlation_bulk_secretory.pdf", width = 7, height = 8)
# save_pheatmap_pdf(pp, "Correlation_bulk_secretory_soupX.pdf", width = 7, height = 8)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
p.mat <- cor.mtest(log_bulk_sample)
corrplot(ct, method = "number", type="upper")
corrplot(ct, method="color", col=col(200),  #cl.lim = c(0.9, 1),
         type="upper", order="hclust", 
         addCoef.col = "black", # Ajout du coefficient de corrélation
         tl.col="black", tl.srt=45, #Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # Cacher les coefficients de corrélation sur la diagonale
         diag=FALSE 
)

## Remove the genes with zeros values despite the count addition
zeros = which(rowSums(bulk_samples)<2)
bulk_samples <- bulk_samples[-zeros,]


## Create a dge object, normalize using TMM
dge <- DGEList(round(bulk_samples), group = serie$type) 
keep <- rowSums(bulk_samples>1) > 5 # 12071, no additional filtering required actually
table(keep)
dge <- dge[keep,]
dge <- calcNormFactors(dge)  # Normalize data

## Create the design matrix (blocking on donor ids)
Group <- factor(serie$type)
design <- model.matrix(~ Group) 
dge <- estimateDisp(dge, design = design) 
fit <- glmFit(dge, design = design) 
colnames(fit$coefficients)

# Differential expression testing
lrt <- glmLRT(fit, coef=2) 
tt <- topTags(lrt, n = Inf) 

# Visual output of the differential expression QCs qnd results
plotBCV(dge) 
hist(tt$table$PValue, 50) 
hist(tt$table$FDR, 50) 
hist(tt$table$logFC)

tt$table$toplot <- rep(0,nrow(tt$table))
tt$table$toplot[tt$table$FDR<0.05 & abs(tt$table$logFC) > 4 & abs(tt$table$logCPM) > 3 ] <- 1
table(tt$table$toplot)
tt$table$Gene = rownames(tt$table)

gg <- ggplot(data=tt$table, aes(x=logCPM,
                                y=logFC)) +  
  geom_point(color="grey", size=1) +
  geom_point(data=subset(tt$table,toplot==1),colour = "red", size=2) +
  theme_bw() +
  geom_text_repel(data=subset(tt$table,toplot==1) ,aes(label=Gene), size=3, 
                  show.legend=FALSE, force=2, colour="blue")
gg

## Add a score (logFC * p-value) to the top differentially expressed genes table
topTable <- tt$table
topTable$score <- topTable$logFC * -log10(topTable$FDR)
topTable$score[topTable$score == -Inf] = 0
topTable$score[topTable$score == Inf] = max(topTable$score) + 1

write.table(topTable, file = "DA_bulk_secretory_NvsB_final.tsv",
            sep = "\t", quote = F, row.names = F, col.names = T)

# write.table(topTable, file = "DA_bulk_secretory_NvsB_soupX_corrected.csv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)


##########################################################################################
#              Differential expression testing : Multiciliated cells                     #
##########################################################################################

## Filter out Ribo genes and MT genes
RB_genes = read.delim("/data/deprez_data/HCA/PeerLab_analysis/RB_genes", header = F, stringsAsFactors = F)
sbt_count = countTable[!(rownames(countTable) %in% 
                           c(grep("^MT-", rownames(countTable), value = T), RB_genes[,"V1"])),]


# Estimate the number of concerned cells
table(metadata$CellType.Corrected)
table(metadata$CellType.Corrected, metadata$Sample)
metadata$annotation <- metadata$CellType.Corrected

## Subset dataset ------------------------------------------------------------------------

# Remove sample with very different cell type composition and thus very different gene expression background
sbt_count <- sbt_count[, rownames(metadata[metadata$annotation %in% c("Multiciliated", "Multiciliated N") &
                                             !metadata$Sample %in% c("D339_Brus_Dis1", "D363_Brus_Dis1"), ])]
sbt_meta <- metadata[rownames(metadata[metadata$annotation %in% c("Multiciliated", "Multiciliated N") &
                                         !metadata$Sample %in% c("D339_Brus_Dis1", "D363_Brus_Dis1"), ]), ]

table(sbt_meta$annotation, sbt_meta$Sample)
table(sbt_meta$annotation)

sbt_meta$bulk_id <- ""
sbt_meta[sbt_meta$annotation == "Multiciliated", "bulk_id"] <- paste0("bulk_MCCs_", sample(
  rep(1:24, 340), sum(sbt_meta$annotation == "Multiciliated"), replace=FALSE))
sbt_meta[sbt_meta$annotation == "Multiciliated N", "bulk_id"] <- paste0("bulk_MCCsN_", sample(
  rep(1:5, 340), sum(sbt_meta$annotation == "Multiciliated N"), replace=FALSE))

# Create bulk samples
bulk_samples <- matrix(data = 0, nrow = nrow(sbt_count),
                       ncol = length(unique(sbt_meta$bulk_id)))
colnames(bulk_samples) <- unique(sbt_meta$bulk_id)

for (sample in unique(sbt_meta$bulk_id)){
  print(sample)
  ix <- rownames(sbt_meta[sbt_meta$bulk_id == sample,])
  print(length(ix))
  bulk_samples[, sample] <- rowSums(sbt_count[, ix])
}
rownames(bulk_samples) <- rownames(sbt_count)


## Prepare summary file - design for statistical testing
serie <- data.frame(Sample=colnames(bulk_samples),
                    type = unlist(strsplit(colnames(bulk_samples), "_"))[seq(2, ncol(bulk_samples)*3, 3)])

## Check for similar library size
serie$Counts <- colSums(bulk_samples, na.rm = T)
gg <-ggplot(serie, aes(Sample, Counts)) + 
  geom_bar(stat = "identity", aes(fill=type)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
gg
# pdf("library_size_bulk_multiciliated.pdf", width=7, height=4, useDingbats=FALSE)
# # pdf("library_size_bulk_multiciliated_soupX.pdf", width=7, height=4, useDingbats=FALSE)
# dev.off()


## Correlation between cluster
mycols=brewer.pal(9,"Spectral")
names(mycols) <- unique(serie$type)
annotation_col = data.frame(Type = serie$type,
                            Position = serie$Sample)
rownames(annotation_col) = serie$Sample

ann_colors = list(Type=mycols)

log_bulk_sample <- log(bulk_samples + 1)
ct <- cor(log_bulk_sample)
pp<- pheatmap(ct, show_rownames = F, show_colnames = T, 
              annotation_col = annotation_col, annotation_colors = ann_colors)
# save_pheatmap_pdf(pp, "Correlation_bulk_multiciliated.pdf", width = 7, height = 8)
# save_pheatmap_pdf(pp, "Correlation_bulk_multiciliated_soupX.pdf", width = 7, height = 8)

## Remove the zeros here
zeros = which(rowSums(bulk_samples)<2)
bulk_samples <- bulk_samples[-zeros,]

## Create a dge object, normalize using TMM
dge <- DGEList(round(bulk_samples), group = serie$type) 
keep <- rowSums(bulk_samples>1) > 5 # 12071, no additional filtering required actually
table(keep)
dge <- dge[keep,]
dge <- calcNormFactors(dge) 

## Create the design matrix (blocking on donor ids)
Group <- factor(serie$type)
design <- model.matrix(~ Group) 
#colnames(design)
dge <- estimateDisp(dge, design = design) 
fit <- glmFit(dge, design = design) 
colnames(fit$coefficients)

## Differential expression testing
lrt <- glmLRT(fit, coef=2) 
tt <- topTags(lrt, n = Inf) 

# Check differential expression QCs qnd results
plotBCV(dge) 
hist(tt$table$PValue, 50) 
hist(tt$table$FDR, 50) 
hist(tt$table$logFC)

tt$table$toplot <- rep(0,nrow(tt$table))
tt$table$toplot[tt$table$FDR<0.05 & abs(tt$table$logFC) >3 & abs(tt$table$logCPM) > 3 ] <- 1
table(tt$table$toplot)
tt$table$Gene = rownames(tt$table)

ggplot(data=tt$table, aes(x=logCPM,
                          y=logFC)) +  
  geom_point(color="grey", size=1) +
  geom_point(data=subset(tt$table,toplot==1),colour = "red", size=2) +
  theme_bw() +
  geom_text_repel(data=subset(tt$table,toplot==1) ,aes(label=Gene), size=3, 
                  show.legend=FALSE, force=2, colour="blue")

topTable <- tt$table
topTable$score <- topTable$logFC * -log10(topTable$FDR)


# write.table(topTable, file = "DA_bulk_multiciliated_final.tsv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)
# write.table(topTable, file = "DA_bulk_multiciliated_soupX_corrected.csv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)



##########################################################################################
#                 Differential expression testing : Suprabasal cells                     #
##########################################################################################

## Filter out Ribo genes and MT genes
RB_genes = read.delim("/data/deprez_data/HCA/PeerLab_analysis/RB_genes", header = F, stringsAsFactors = F)
sbt_count = countTable[!(rownames(countTable) %in% 
                           c(grep("^MT-", rownames(countTable), value = T), RB_genes[,"V1"])),]

metadata$annotation <- metadata$CellType.Corrected

## Subset dataset
sbt_count <- sbt_count[, rownames(metadata[metadata$annotation %in% c("Suprabasal", "Suprabasal N") &
                                             !metadata$Sample %in% c("D322_Biop_Nas1", "D339_Biop_Int1"), ])]
sbt_meta <- metadata[rownames(metadata[metadata$annotation %in% c("Suprabasal", "Suprabasal N") &
                                         !metadata$Sample %in% c("D322_Biop_Nas1", "D339_Biop_Int1"), ]), ]

table(sbt_meta$annotation, sbt_meta$Sample)
table(sbt_meta$annotation)

sbt_meta$bulk_id <- ""
sbt_meta[sbt_meta$annotation == "Suprabasal", "bulk_id"] <- paste0("bulk_SB_", sample(
  rep(1:11, 700), sum(sbt_meta$annotation == "Suprabasal"), replace=FALSE)) #11

# soupX
# sbt_meta[sbt_meta$annotation == "Suprabasal", "bulk_id"] <- paste0("bulk_SB_", sample(
#   rep(1:14, 570), sum(sbt_meta$annotation == "Suprabasal"), replace=FALSE))
sbt_meta[sbt_meta$annotation == "Suprabasal N", "bulk_id"] <- paste0("bulk_SBN_", sample(
  rep(1:9, 470), sum(sbt_meta$annotation == "Suprabasal N"), replace=FALSE))


bulk_samples <- matrix(data = 0, nrow = nrow(sbt_count),
                       ncol = length(unique(sbt_meta$bulk_id)))
colnames(bulk_samples) <- unique(sbt_meta$bulk_id)

for (sample in unique(sbt_meta$bulk_id)){
  print(sample)
  ix <- rownames(sbt_meta[sbt_meta$bulk_id == sample,])
  print(length(ix))
  bulk_samples[, sample] <- rowSums(sbt_count[, ix])
}
rownames(bulk_samples) <- rownames(sbt_count)

## Prepare summary file - design for statistical testing
serie <- data.frame(Sample=colnames(bulk_samples),
                    type = unlist(strsplit(colnames(bulk_samples), "_"))[seq(2, ncol(bulk_samples)*3, 3)])

## Library size
serie$Counts <- colSums(bulk_samples, na.rm = T)
gg <-ggplot(serie, aes(Sample, Counts)) + 
  geom_bar(stat = "identity", aes(fill=type)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
gg
# pdf("library_size_bulk_suprabasal.pdf", width=7, height=4, useDingbats=FALSE)
# # pdf("library_size_bulk_suprabasal_soupX.pdf", width=7, height=4, useDingbats=FALSE)
# dev.off()

## Correlation between cluster
mycols=brewer.pal(9,"Spectral")
names(mycols) <- unique(serie$type)
annotation_col = data.frame(Type = serie$type,
                            Position = serie$Sample)
rownames(annotation_col) = serie$Sample
ann_colors = list(Type=mycols)

log_bulk_sample <- log(bulk_samples + 1)
ct <- cor(log_bulk_sample)
pp<- pheatmap(ct, show_rownames = F, show_colnames = T, 
              annotation_col = annotation_col, annotation_colors = ann_colors)
# save_pheatmap_pdf(pp, "Correlation_bulk_multiciliated.pdf", width = 7, height = 8)
# save_pheatmap_pdf(pp, "Correlation_bulk_multiciliated_soupX.pdf", width = 7, height = 8)

## Remove the zeros here
zeros = which(rowSums(bulk_samples)<2)
bulk_samples <- bulk_samples[-zeros,]

## Create a dge object, normalize using TMM
dge <- DGEList(round(bulk_samples), group = serie$type) 
keep <- rowSums(bulk_samples>1) > 5 # 12071, no additional filtering required actually
table(keep)
dge <- dge[keep,]
dge <- calcNormFactors(dge) 

## Create the design matrix (blocking on donor ids)
Group <- factor(serie$type)
design <- model.matrix(~ Group) 
dge <- estimateDisp(dge, design = design) 
fit <- glmFit(dge, design = design) 
colnames(fit$coefficients)


## Differential expression testing
lrt <- glmLRT(fit, coef=2) 
tt <- topTags(lrt, n = Inf) 

# Differential expression QCs qnd results
plotBCV(dge) 
hist(tt$table$PValue, 50) 
hist(tt$table$FDR, 50) 
hist(tt$table$logFC)

tt$table$toplot <- rep(0,nrow(tt$table))
tt$table$toplot[tt$table$FDR<0.05 & abs(tt$table$logFC) > 4 & abs(tt$table$logCPM) > 4 ] <- 1
table(tt$table$toplot)
tt$table$Gene = rownames(tt$table)

gg <- ggplot(data=tt$table, aes(x=logCPM,
                                y=logFC)) +  
  geom_point(color="grey", size=1) +
  geom_point(data=subset(tt$table,toplot==1),colour = "red", size=2) +
  theme_bw() +
  geom_text_repel(data=subset(tt$table,toplot==1) ,aes(label=Gene), size=3, 
                  show.legend=FALSE, force=2, colour="blue")
gg + ggtitle("MA-plot: Proximal biopsies - bulk")

topTable <- tt$table
topTable$score <- topTable$logFC * -log10(topTable$FDR)


# write.table(topTable, file = "DA_bulk_suprabasal.tsv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)

# write.table(topTable, file = "DA_bulk_suprabasal_soupX_corrected.csv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)



##########################################################################################
#                              Gene Set Enrichment Analysis                              #
##########################################################################################

data(kegg.gs)

pathways_kegg <- gmtPathways("/home/deprez/Misc/Gsea/c2.cp.kegg.v6.2.symbols.gmt")
pathways_cp <- gmtPathways("/home/deprez/Misc/Gsea/c2.cp.v6.0.symbols.gmt") # canonical pathways
pathways_reactome <- gmtPathways("/home/deprez/Misc/Gsea/c2.cp.reactome.v6.2.symbols.gmt")
pathways_GO_bp <- gmtPathways("/home/deprez/Misc/Gsea/c5.bp.v6.2.symbols.gmt")
pathways_GO_cc <-gmtPathways("/home/deprez/Misc/Gsea/c5.cc.v6.2.symbols.gmt")
pathways_GO_mf <- gmtPathways("/home/deprez/Misc/Gsea/c5.mf.v6.2.symbols.gmt")

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

topTable <- read.table("DA_bulk_secretory.csv", 
                       sep = "\t", header = T)
data <- topTable$logFC
names(data) <- topTable$Gene
data <- sort(data, decreasing = T)
data = data[0:1000]
data[1:10]
fgseaRes_GO_bp <- fgsea(pathways = pathways_GO_bp, stats = data, minSize=15, maxSize=500, nperm=100000)


df_GO_bp <- data.frame(fgseaRes_GO_bp)
df_GO_bp$leadingEdge <- ""
df_GO_bp$Fraction <- 0
df_GO_bp$LE_size <- 0

for(i in 1:nrow(fgseaRes_GO_bp)){
  df_GO_bp$leadingEdge[i] <- paste(fgseaRes_GO_bp$leadingEdge[[i]], sep = ",", collapse = ",")
  df_GO_bp$LE_size[i] <- length(unlist(fgseaRes_GO_bp$leadingEdge[i]))
  df_GO_bp$Fraction[i] <- round((df_GO_bp$LE_size[i] / fgseaRes_GO_bp$size[i]) * 100)
}


# write.table(df_GO_bp, file = "./GSEA/secretory_NASAL_fgseaRes_GO_bp.csv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)
# write.table(df_GO_bp, file = "./GSEA/secretory_BRONCHIAL_fgseaRes_GO_bp.csv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

topTable <- read.table("DA_bulk_suprabasal.csv", 
                       sep = "\t", header = T)
data <- topTable$logFC
names(data) <- topTable$Gene
data <- sort(data, decreasing = T)
data = data[0:500]
data[1:10]
fgseaRes_GO_bp <- fgsea(pathways = pathways_GO_bp, stats = data, minSize=15, maxSize=500, nperm=100000)

df_GO_bp <- data.frame(fgseaRes_GO_bp)
df_GO_bp$leadingEdge <- ""
df_GO_bp$Fraction <- 0
df_GO_bp$LE_size <- 0

for(i in 1:nrow(fgseaRes_GO_bp)){
  df_GO_bp$leadingEdge[i] <- paste(fgseaRes_GO_bp$leadingEdge[[i]], sep = ",", collapse = ",")
  df_GO_bp$LE_size[i] <- length(unlist(fgseaRes_GO_bp$leadingEdge[i]))
  df_GO_bp$Fraction[i] <- round((df_GO_bp$LE_size[i] / fgseaRes_GO_bp$size[i]) * 100)
}

# write.table(df_GO_bp, file = "./GSEA/suprabasal_NASAL_fgseaRes_GO_bp.csv",
#             sep = "\t", quote = F, row.names = F, col.names = T)
# write.table(df_GO_bp, file = "./GSEA/suprabasal_BRONCHIAL_fgseaRes_GO_bp.csv",
#             sep = "\t", quote = F, row.names = F, col.names = T)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

topTable <- read.table("DA_bulk_multiciliated.csv", 
                       sep = "\t", header = T)
data <- topTable$logFC
names(data) <- topTable$Gene
data <- sort(data, decreasing = T)
data = data[0:1000]
data[1:10]
fgseaRes_GO_bp <- fgsea(pathways = pathways_GO_bp, stats = data, minSize=15, maxSize=500, nperm=100000)

df_GO_bp <- data.frame(fgseaRes_GO_bp)
df_GO_bp$leadingEdge <- ""
for(i in 1:nrow(fgseaRes_GO_bp)){
  df_GO_bp$leadingEdge[i] <- paste(fgseaRes_GO_bp$leadingEdge[[i]], sep = ",", collapse = ",")
}
head(df_GO_bp$pathway)

# write.table(df_GO_bp, file = "./GSEA/multiciliated_NASAL_fgseaRes_GO_bp.csv",
#             sep = "\t", quote = F, row.names = F, col.names = T)
# write.table(df_GO_bp, file = "./GSEA/multiciliated_BRONCHIAL_fgseaRes_GO_bp.csv",
#             sep = "\t", quote = F, row.names = F, col.names = T)


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# fgseaRes_GO_cc <- fgsea(pathways = pathways_GO_cc, stats = data, minSize=15, maxSize=500, nperm=100000)
# fgseaRes_GO_mf <- fgsea(pathways = pathways_GO_mf, stats = data, minSize=15, maxSize=500, nperm=100000)
# fgseaRes_react <- fgsea(pathways = pathways_reactome, stats = data, minSize=15, maxSize=500, nperm=100000)
# fgseaRes_kegg <- fgsea(pathways = pathways_kegg, stats = data, minSize=15, maxSize=500, nperm=100000)

# df_react <- data.frame(fgseaRes_react)
# df_react$leadingEdge <- ""
# for(i in 1:nrow(fgseaRes_react)){
#   df_react$leadingEdge[i] <- paste(fgseaRes_react$leadingEdge[[i]], sep = ",", collapse = ",")
# }
# 
# df_kegg <- data.frame(fgseaRes_kegg)
# df_kegg$leadingEdge <- ""
# for(i in 1:nrow(fgseaRes_kegg)){
#   df_kegg$leadingEdge[i] <- paste(fgseaRes_kegg$leadingEdge[[i]], sep = ",", collapse = ",")
# }
# 
# df_GO_cc <- data.frame(fgseaRes_GO_cc)
# df_GO_cc$leadingEdge <- ""
# for(i in 1:nrow(fgseaRes_GO_cc)){
#   df_GO_cc$leadingEdge[i] <- paste(fgseaRes_GO_cc$leadingEdge[[i]], sep = ",", collapse = ",")
# }
# 
# df_GO_mf <- data.frame(fgseaRes_GO_mf)
# df_GO_mf$leadingEdge <- ""
# for(i in 1:nrow(fgseaRes_GO_mf)){
#   df_GO_mf$leadingEdge[i] <- paste(fgseaRes_GO_mf$leadingEdge[[i]], sep = ",", collapse = ",")
# }

# write.table(df_react, file = "secretoryN_fgseaRes_react.csv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)
# write.table(df_kegg, file = "fgseaRes_kegg.csv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)
# write.table(df_GO_cc, file = "secretoryN_fgseaRes_GO_cc.csv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)
# write.table(df_GO_mf, file = "secretoryN_fgseaRes_GO_mf.csv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)

