import scanpy as sc


def preprocess(adata,min_genes=400,batch=False,mit_cutoff=0.25,max_counts=40000):
    sc.pp.filter_cells(adata,min_genes=min_genes)
    sc.pp.filter_cells(adata,max_counts=max_counts)
    #adata = adata[adata.obs['pct_mito'] < mit_cutoff]
    adata.raw = adata.copy()
    sc.pp.highly_variable_genes(adata,flavor='seurat_v3',n_top_genes=3000, subset=True)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata, copy=False)
    adata.layers['log1p'] = adata.X.copy()
    sc.pp.scale(adata,copy=False)
    sc.pp.pca(adata,n_comps=30)
    if batch:
        sc.external.pp.harmony_integrate(adata,key='sample',basis='X_pca',adjusted_basis='X_pca_harmony')
    return adata