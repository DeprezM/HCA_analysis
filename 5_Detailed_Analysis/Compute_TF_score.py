#!usr/bin/python

# Import useful modules

import numpy as np
import pandas as pd
import scanpy as sc
import os
#import igraph
import matplotlib.pyplot as plt
import seaborn
import logging as logg
from scipy.sparse import issparse
import sys


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=90)


########################################################################################################
########################################################################################################
def score_genes_to_list(
        adata,
        gene_list,
        ctrl_size=50,
        n_bins=25,
        random_state=0):  # we use the scikit-learn convention of calling the seed "random_state"
    
    """Score a set of genes [Satija15]_.
    Imported and modify from scanpy _ function.
    The score is the average expression of a set of genes subtracted with the
    average expression of a reference set of genes. The reference set is
    randomly sampled from the `gene_pool` for each binned expression value.
    This reproduces the approach in Seurat [Satija15]_ and has been implemented
    for Scanpy by Davide Cittaro.
    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix.
    gene_list : iterable
        The list of gene names used for score calculation.
    ctrl_size : `int`, optional (default: 50)
        Number of reference genes to be sampled. If `len(gene_list)` is not too
        low, you can set `ctrl_size=len(gene_list)`.
    n_bins : `int`, optional (default: 25)
        Number of expression level bins for sampling.
    random_state : `int`, optional (default: 0)
        The random seed for sampling.
        
    Returns
    -------
    A list with the computed score.
    Examples
    --------
    See this `notebook <https://github.com/theislab/scanpy_usage/tree/master/180209_cell_cycle>`__.
    """
    
    if random_state:
        np.random.seed(random_state)

    gene_list_in_var = []
    var_names = adata.var_names
    for gene in gene_list:
        if gene in var_names:
            gene_list_in_var.append(gene)
        else:
            logg.warning(f'gene: {gene} is not in adata.var_names and will be ignored')
    gene_list = set(gene_list_in_var[:])

    
    gene_pool = list(var_names)

    # Trying here to match the Seurat approach in scoring cells.
    # Basically we need to compare genes against random genes in a matched
    # interval of expression.

    # TODO: this densifies the whole data matrix for `gene_pool`
    if issparse(adata.X):
        obs_avg = pd.Series(
            np.nanmean(
                adata[:, gene_pool].X.toarray(), axis=0), index=gene_pool)  # average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(adata[:, gene_pool].X, axis=0), index=gene_pool)  # average expression of genes

    obs_avg = obs_avg[np.isfinite(obs_avg)] # Sometimes (and I don't know how) missing data may be there, with nansfor

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method='min') // n_items
    control_genes = set()

    # now pick `ctrl_size` genes from every cut
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = np.array(obs_cut[obs_cut == cut].index)
        np.random.shuffle(r_genes)
        control_genes.update(set(r_genes[:ctrl_size]))  # uses full r_genes if ctrl_size > len(r_genes)

    # To index, we need a list - indexing implies an order.
    control_genes = list(control_genes - gene_list)
    gene_list = list(gene_list)


    X_list = adata[:, gene_list].X
    if issparse(X_list): X_list = X_list.toarray()
    X_control = adata[:, control_genes].X
    if issparse(X_control): X_control = X_control.toarray()
    X_control = np.nanmean(X_control, axis=1)


    score = np.nanmean(X_list, axis=1) - X_control
    
    del X_list, X_control, gene_list_in_var, var_names, gene_list, ctrl_size,
    del gene_pool, obs_avg, n_items, obs_cut, control_genes, r_genes
    
    return np.array(score).ravel().tolist()

########################################################################################################
########################################################################################################


## Load datasets
adata = sc.read_h5ad('/Data/Annotated_dataset_v1.h5ad')
tf_table = pd.read_csv('/Data/table_top50_TF.tsv', sep = "\t", 
                      header = 0, index_col= 0)

# Create 0 filled dataframe
#result_array = np.empty(shape=(len(tf_table.columns.tolist()),len(adata.obs.index.tolist())))
#tf_score = pd.DataFrame(a, index = tf_table.columns.tolist(), 
#                        columns = adata.obs.index.tolist())

if int(sys.argv[1]) == 0:
  file_object = open('/data/deprez_data/HCA/Analysis/FullDataset_v4/score_table_top50_TF.tsv', 'w') 
else :
  file_object = open('/data/deprez_data/HCA/Analysis/FullDataset_v4/score_table_top50_TF.tsv', 'a') 

  
  
if int(sys.argv[1]) == 0:
    file_object.write('TF\t' + '\t'.join(adata.obs.index.tolist()) + '\n')

# Fill in the table
for tf in tf_table.columns.tolist()[int(sys.argv[1]):int(sys.argv[1])+5] :
    print(tf)

    result = score_genes_to_list(adata, tf_table.loc[:, tf].tolist())
    file_object.write(tf + '\t' + '\t'.join(map(str,result)) + '\n')


file_object.close()



