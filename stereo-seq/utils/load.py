import os
import scanpy as sc
import anndata as ad

def merge_dfs(samples,f="filtered_feature_bc_matrix.h5",joint_dir="data/aggregated",fname="aggregated.h5ad"):
    adata_arr = []
    for sample in samples:
        f_f = os.path.join(DATA_DIR,sample,"outs",f)
        assert os.path.exists(f_f)
        adata = sc.read_10x_h5(f_f)
        adata.var_names_make_unique()
        adata_arr.append(adata)
    
    merged_ann = ad.concat(adata_arr,axis=0,label='sample', index_unique='_', keys =samples)
    out_f = os.path.join(joint_dir,fname)
    merged_ann.write_h5ad(out_f)
    return merged_ann