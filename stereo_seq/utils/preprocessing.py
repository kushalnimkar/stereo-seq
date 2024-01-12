import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


import os
import re


def pipeline(adata,batch=True,mit_cutoff=0.25,min_genes=400,max_counts=40000,output_f="",figure_dir=""):
    if batch:
        pre="batch_correct"
        use_rep='X_pca_harmony'
    else:
        pre="no_batch_correct"
        use_rep='X_pca'
        mit_cutoff=1.0

    precluster_QC(adata)
    adata = preprocess(adata,batch=batch,mit_cutoff=mit_cutoff,min_genes=min_genes,max_counts=max_counts)

    cluster_and_umap(adata,use_rep=use_rep)   
    make_QC_plots(adata,figure_dir,pre=pre)

    adata.write_h5ad(output_f)
    return adata

def precluster_QC(adata):
    pat = re.compile(r'mt')
    mito_genes = [i for i in adata.var.index if pat.match(i)]
    adata.obs['mRNA_count'] = np.sum(adata.X,axis=1)
    adata.obs['n_mito'] = ( np.array(np.sum(adata[:,mito_genes].X,axis=1)).reshape(-1) )

    adata.obs['pct_mito'] = adata.obs['n_mito']/adata.obs['mRNA_count']
    adata.obs['n_genes'] = np.array(np.sum(adata.X>0,axis=1)).reshape(-1)

    return mito_genes

def cluster_and_umap(adata,use_rep, cluster=True,umap=True):
    if cluster:
        sc.pp.neighbors(adata,n_neighbors=50,use_rep=use_rep)
        sc.tl.leiden(adata)
    if umap:
        sc.tl.umap(adata)
    return adata

def preprocess(adata,min_genes=400,batch=False,mit_cutoff=0.25,fixed_genes=None,max_counts=40000):
    adata=adata.copy()
    sc.pp.filter_cells(adata,min_genes=min_genes)
    sc.pp.filter_cells(adata,max_counts=max_counts)

    if 'pct_mito' in adata.obs.columns:
        adata = adata[adata.obs['pct_mito'] < mit_cutoff]
    else:
        print("Note: Not filtering for cells with too many mitochondrial genes")

    adata.raw = adata.copy()
    if fixed_genes is not None:
        adata._inplace_subset_var(fixed_genes)
    else:
        sc.pp.highly_variable_genes(adata,flavor='seurat_v3',n_top_genes=3000, subset=False)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata, copy=False)
    adata.layers['log1p'] = adata.X.copy()
    sc.pp.scale(adata,copy=False)
    sc.pp.pca(adata,n_comps=30)
    if batch:
        sc.external.pp.harmony_integrate(adata,key='sample',basis='X_pca',adjusted_basis='X_pca_harmony')
    return adata


def make_QC_plots(adata,figure_dir,pre):
    with plt.rc_context():
        
        sc.pl.umap(adata,color='sample',show=False)
        plt.savefig(os.path.join(figure_dir,f"{pre}_samples.png"))
        plt.show()

        sc.pl.umap(adata,color='pct_mito',vmin=0,vmax=0.2, show=False)
        plt.savefig(os.path.join(figure_dir,f"{pre}_mito.png"))
        plt.show()

        sc.pl.umap(adata,color='pct_mito',vmin=0,vmax=1, show=False)
        plt.savefig(os.path.join(figure_dir,f"{pre}_mito_01.png"))
        plt.show()

        plt.hist(adata.obs['pct_mito'],cumulative=True,density=True,bins=1000)
        plt.xlabel("pct_mito")
        plt.ylabel("Fraction of cells < pct_mito")
        plt.savefig(os.path.join(figure_dir,f"{pre}_pct_mito.png"),dpi=200,bbox_inches='tight')
        plt.show()


def split_adata(adata,groups,output_dir,other_group='Non-neuronal',cluster_col='leiden',class_col='class',output_pre='NR_'):
    df_l = []
    g_lst=[]
    for group in groups:
        g_lst.append(group)
        info = groups[group]
        assert len(info['marker_genes'])== len(info['z-scores'])
        cutoffs = np.array(adata[:,info['marker_genes']].X)> np.array(info['z-scores'])

        adata.obs[group]=np.any(cutoffs,axis=1)
        agg = adata.obs.groupby(cluster_col)[group].mean().reset_index(name=group)
        agg.set_index(cluster_col,inplace=True)
        plt.hist(agg[group])
        plt.ylabel('Number of clusters')
        plt.xlabel('Prop of cells above cutoff zscore for at least 1 {} marker gene'.format(group))
        plt.show()

        potential_clusters = agg[agg[group] >info['pct_in']].index.values
        adata.obs[group] = pd.Categorical(adata.obs[cluster_col].isin(potential_clusters))

        sc.pl.umap(adata,color=[cluster_col]+info['marker_genes'] + [group],vmin=0,vmax=2)



        df_l.append(agg)
    
    bool_t_ident = adata.obs[g_lst].values
    assert not np.any(np.sum(bool_t_ident,axis=1) >= 2),"One cluster is assigned to two classes"
    
    idx_identities = ( np.sum(bool_t_ident,axis=1) !=0 ) * (bool_t_ident.argmax(axis=1)+1) #makes 0 the none type, 1 excitatory, 2 inhibitory...
    mapping = dict( zip(range(len(g_lst)+1), [other_group] + g_lst ) )
    adata.obs[class_col] = pd.Categorical( pd.Series(idx_identities,index=adata.obs.index).apply(lambda x: mapping[x]) )

    for cat in adata.obs[class_col].cat.categories:
        adata_c = adata[adata.obs[class_col]==cat].copy()
        output_f = os.path.join(output_dir,"{}{}.h5ad".format(output_pre,cat))
        adata_c.X = adata_c.raw.X
        adata_c.write_h5ad(output_f)
    
    return pd.concat(df_l,axis=1) 
    

    

    

