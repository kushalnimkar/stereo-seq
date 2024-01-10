import scanpy as sc
import matplotlib.pyplot as plt
import re
import numpy as np


import os



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

def preprocess(adata,min_genes=400,batch=False,mit_cutoff=0.25,max_counts=40000):
    sc.pp.filter_cells(adata,min_genes=min_genes)
    sc.pp.filter_cells(adata,max_counts=max_counts)
    adata = adata[adata.obs['pct_mito'] < mit_cutoff]
    print(adata)
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

def OLD_make_plots(adata,figure_dir):
    plt.scatter(adata.obs['n_genes'],adata.obs['pct_mito'],s=0.1)
    plt.xlabel("Number of genes with non-zero expression",fontsize=14)
    plt.ylabel("Fraction of total mRNA that is mitochondrial gene",fontsize=14)
    plt.show()
    plt.hist(adata.obs['pct_mito'],cumulative=True,density=True,bins=1000)
    plt.xlabel("pct_mitochondrial")
    plt.ylabel("Fraction of all cells")
    plt.title("Cumulative Distribution Function")
    plt.show()

    plt.hist(adata.obs['n_genes'],cumulative=True,density=True,bins=1000)
    plt.xlabel('n_genes')
    plt.ylabel("Fraction of all cells")
    plt.title("Cumulative Distribution Function")
    plt.show()

    

    

    

