# Analysis of the healthy human airways at single-cell level.
The code is structured in 6 steps, each with specific inputs and outputs (available in the data section).

### 1. Primary data analysis
Example code of the individual exploratory analysis performed on each of the 35 samples composing the atlas, with use of the Seurat v3 R package. From this initial analysis, repeated on the 35 samples, we obtained a rough estimate of the cell type composition of the atlas (to be refined later in the analysis of the complete dataset), and a list of robust marker genes found expressed in each cell types identified across all 35 samples.

**Output from this step :** `Robust_marker_genes.tsv`
**Example** `D344_Biop_Int1`

### 2. Pre-processing of the data
Pre-processed of the complete dataset in parallel in 3 different scripts: 

   - Consensus cells and genes filtering across the 35 samples, normalisation of all the cell count to 10000 UMIs and merging of all the samples in a single dataset producing a large and preliminary-processed count table. [`Pre-processing_preliminary_dataset.ipynb`] 

   - Iterative preliminary analysis of the dataset includes progressive cell filtering of small clusters composed of 'low quality'/'peculiar' cells (high mitochondrial cluster-cells, ...). [`Preliminary analysis_v1 ...v4`]
    
   - Consensus cells- and genes-filtering across the 35 samples (without normalisation) and merging of all the datasets in a single one to produce a large and unique raw count table. [`PreProcessing_raw_dataset.ipynb`]

   - Identification/Inference of doublet cells across all the 35 samples independently, and further analysis of the dataset to estimate the proportion of inferred doublet cells in the resulting clusters and corresponding cell filtering. [`PreProcessing_doublets.ipynb`; `Pre-analysis_doublets.ipynb`; `Preliminary_analysis_doublet_metadata.Rmd`]   
    
   - Pre-processing of the background in gene expression across all samples to produce a 'background free' raw count table. [`Pre-processing_gene background_dataset.Rmd`; `Preliminary_background_analysis`]
   
**Input files for this step**
All the 10x output files from the 35 samples (available for download on GSE)
RB_genes (list of the filtered out ribosomal genes)
   
**Output files from this step:**
Preliminary Analysis datasets :
  - `PreProcessed_preliminary_dataset`
  - `Preliminary_analysis_v1...v4`

PreProcessed raw count Table:
  - `PreProcessed_raw_dataset`
  
Doublet analysis:
  - `metadata_doublet`
  
Background free datasets :
  - `background_features`
  - `background_metadata`
  - `SoupX_raw_dataset`
  - `SoupX_strained_dataset` (Advance SoupX correction, not used in the following analysis)
  
### 3. Data Normalization and integration (batch correction)

To appropriately normalize the complete dataset, it is again filtered of the 'low-quality cells' and then normalised using scran method (Lun & Haghverdi) on both raw counts and soupX corrected counts. Lastly, the normalised counts data are integrated to produce a batch-free PCA matrix that will be used in the following analysis.

**Output files**
  - `scranNorm_dataset`
  - `fastMNN_PCA`
  - `SoupX_norm_dataset`

**Integration process** 
The data integration process will progressively map one dataset onto another in the following order :
Intermediate samples, Tracheal (Proximal) samples, Distal samples and Nasal samples. The order was defined based on the results of the preliminary analysis, which established the relative homogeneity of the samples based on their sampling location. For all the samples from the 'same' initial sampling site (level), they are aggregated from the larger to the smaller datasets (the ones with more cells first).


### 4. Analysis

The dataset can now be fully analysed including umap embedding (computed on the integrated PCA), clustering of the cells, marker genes identification and specific sub-clustering of each 'key' cell-cluster. The many sub-clustering steps were done to improve the precision of the final cell labelling [`Annotated_dataset_metadata`].

**Output files**
  -`Annotated_dataset_v1` (v1 because the cell types names will be progressively updated as our understanding improves)
  -`Focus_XXX_cells`
  -`markers_XXX_cells`

### 5. Detailed_Analysis

This repertory contains the scripts used for the detailed analysis of some cell types. It includes :

   - Differential analysis of the similar cell types identified in both Nasal and Bronchial samples (Secretory, Multiciliated, Suprabasal), followed by Gene Set Enrichment Analysis;
   - Trajectory inference of the epithelial cell types from the Nasal or Bronchial area using PAGA;
   - Inference of the Transcription Factors activity using SCENIC to identify the regulons.
   
**Output files**
  -`DA_bulk_XXX`
  -`XXX_fgseaRes_GO_bp`
  -`anndata_v6_Paga`
  -`score_TF`

### 6. Figures

All the scripts used for the Figures found in the paper Deprez et al. Scripts are labelled by figure numbers. All the necessary files are in the data repository.

**Final AnnData object and metadata : **
  - `Annotated_dataset.h5ad`
  - `Annotated_dataset_metadata.tsv`
  
**Count tables: raw or normalized**
  -`raw_exprMatrix.Rda`
  -`SoupX_raw_dataset.Rda`
  -`scranNorm_dataset.Rda`
  -`SoupX_norm_dataset.Rda`
  
  

