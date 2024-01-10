import os
import scanpy as sc
import anndata as ad
import pandas as pd
import scipy as sci

OUTS_D="outs"

def merge_dfs(samples,f="filtered_feature_bc_matrix.h5",out_f="aggregated.h5ad",format="mtx"):
    """Merges sample datasets
    Args:
        samples (dict): list of sample_names as keys, corresponding directory names as values
        f (str): Data matrix name to load
        out_f (str): Output h5ad that contains all samples concatenated
        format (str): Choice of file format we are trying to load (3 options supported currently...)
    Writes:
        Aggregated dataframe with all datasets chosen at path out_f
    Returns:
        Same aggregated dataframe that was written
    """

    adata_arr = []
    sample_keys =[]

    for sample_name in samples:

        sample_p = samples[sample_name]
        sample_keys.append(sample_name)

        f_f = os.path.join(sample_p,OUTS_D,f)
        assert os.path.exists(f_f),f_f

        if format=="mtx":
            adata = sc.read_10x_mtx(f_f)
        elif format=="h5":
            adata= sc.read_10x_h5(f_f)
        elif format=="custom":
            adata= custom_load(f_f)
        else:
            raise AssertionError("Improper loading format choice given")

        print(adata.var)
        adata.var_names_make_unique()
        adata_arr.append(adata)
    
    merged_ann = ad.concat(adata_arr,axis=0,label='sample', index_unique='_', keys =sample_keys)
    merged_ann.write_h5ad(out_f)
    return merged_ann


def custom_load(data_folder):
    adata = sc.read_mtx(os.path.join(data_folder,'matrix.mtx.gz'))
    print(adata)
    adata_bc=pd.read_csv(os.path.join(data_folder,'barcodes.tsv.gz'),sep='\t',header=None,index_col=0)
    adata_features=pd.read_csv(os.path.join(data_folder,'features.tsv.gz'),sep='\t',header=None)
    adata= adata.T
    adata.obs.index=adata_bc.index
    adata.obs.index.name = 'cell_id'
    adata.var['ensembl_id']= adata_features[0].tolist()  
    adata.var.index= adata_features[1].tolist()
    adata.var.index.name='gene_name'
    return adata