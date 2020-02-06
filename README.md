# Analysis of the human healthy airways at single-cell level.
The code is structured in 6 steps each with specific inputs and outputs (available in the data section).

### 1. Primary data analysis
Example code of the individual exploratory analysis performed on each of the 35 samples composing the atlas, with use of the Seurat v3 R package. From this initial analysis, repeated on the 35 samples, we obtained a rough estimate of the cell type composition of the atlas (to be refined later in the code) and a list of robust marker genes found expressed in the each cell types identified across all 35 samples.

Output from this step : Robust_marker_genes.tsv

### 2. Preprocessing of the data
In this section of the workflow, the data was pre-processed in paralle in 3 different scripts: 

   - Consensus cells and genes filtering across the 35 samples and merge of all the datasets in a single one to produce a large and unique raw count table. Additional cell filtering on the aggregated count matrix to remove doublets and low quality cells identified after an initial analysis of the complete dataset.
    
   - Identification/Inference of doublet cells across all the 35 samples independantly, and further analysis of the dataset to estimate the proportion of inferred doublet cells in the resulting cluster.   
    
   - Pre-processing of the background in gene expression across all samples to produce a 'background free' raw count table.  


