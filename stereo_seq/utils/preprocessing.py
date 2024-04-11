import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scrublet as scr
import anndata as ad

import stereo_seq.utils.exploration as sue

import os
import re


def pipeline(
    adata,batch=True,batch_key='sample',mit_cutoff=0.25,min_genes=100,max_counts=40000,key_added='leiden',
    output_f="",figure_dir="",make_qc=True,detect_doublets=False,override_info=None,inplace=False
    ):
    if batch:
        pre="batch_correct"
        use_rep='X_pca_harmony'
    else:
        pre="no_batch_correct"
        use_rep='X_pca'
        mit_cutoff=1.0

    precluster_QC(adata)
    if not inplace:
        adata = preprocess(adata,batch=batch,mit_cutoff=mit_cutoff,min_genes=min_genes,max_counts=max_counts,inplace=inplace)
    else:
        preprocess(adata,batch=batch,mit_cutoff=mit_cutoff,min_genes=min_genes,max_counts=max_counts,inplace=inplace)

    cluster_and_umap(adata,use_rep=use_rep,key_added=key_added)
    if make_qc:
        make_QC_plots(adata,figure_dir,pre=pre)
    if detect_doublets:
        run_scrublet(adata,rep_col=batch_key,cluster_key=key_added,recalc=True, score_cutoff=None,override_info=override_info)



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

def cluster_and_umap(adata,use_rep, key_added='leiden',cluster=True,umap=True):
    if cluster:
        sc.pp.neighbors(adata,n_neighbors=50,use_rep=use_rep)
        sc.tl.leiden(adata,key_added=key_added)
    if umap:
        sc.tl.umap(adata)
    return adata

def preprocess(adata,min_genes=400,batch=False,batch_key='sample',mit_cutoff=0.25,fixed_genes=None,max_counts=40000,find_hvg=True,do_pca=True,inplace=False):
    if not inplace:
        adata=adata.copy()
    if min_genes is not None:
        sc.pp.filter_cells(adata,min_genes=min_genes)
    if max_counts is not None:
        sc.pp.filter_cells(adata,max_counts=max_counts)
        
    if 'pct_mito' in adata.obs.columns:
        adata._inplace_subset_obs((adata.obs['pct_mito']<mit_cutoff).index)
    else:
        print("Note: Not filtering for cells with too many mitochondrial genes")

    adata.raw = adata.copy()
    if fixed_genes is not None:
        adata._inplace_subset_var(fixed_genes)
    else:
        if find_hvg:
            sc.pp.highly_variable_genes(adata,flavor='seurat_v3',n_top_genes=3000, subset=False)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata, copy=False)
    adata.layers['log1p'] = adata.X.copy()
    sc.pp.scale(adata,copy=False)
    if do_pca:
        sc.pp.pca(adata,n_comps=30)

    if batch:
        sc.external.pp.harmony_integrate(adata,key=batch_key,basis='X_pca',adjusted_basis='X_pca_harmony')

    if not inplace:
        return adata

def run_scrublet(adata,rep_col='sample',cluster_key='leiden',recalc=True, score_cutoff=None,override_info=None):
    #calc doublets
    if recalc:
        adata.obs['doublet_score'] = pd.NA
        adata.obs['predicted_doublet'] = pd.NA

        for rep in adata.obs[rep_col].unique():
            adata_r = adata[adata.obs[rep_col]==rep]
            scrub = scr.Scrublet(adata_r.raw.X)
            doublet_scores, predicted_doublets = scrub.scrub_doublets()
            adata_r.obs['doublet_score'] = pd.Series(doublet_scores,index=adata_r.obs.index)
            if score_cutoff is None:
                adata_r.obs['predicted_doublet'] = pd.Series(predicted_doublets,index=adata_r.obs.index)
            else:
                adata_r.obs['predicted_doublet'] = pd.Series(doublet_scores > score_cutoff,index=adata_r.obs.index)
            sc.pl.umap(adata_r,color='doublet_score')
            plt.show()
            adata.obs.update(adata_r.obs,overwrite=False)
    
    adata.obs['predicted_doublet']= adata.obs['predicted_doublet'].astype(bool)
    adata.obs['doublet_score']= adata.obs['doublet_score'].astype(float)
    mapping = sue.gen_assignment(adata,col_0=cluster_key,col_1='predicted_doublet',assignment_strat="max")
    [print("Potential doublet cluster", x,mapping[x]) for x in mapping if mapping[x]==True]
    adata.obs['is_doublet'] = adata.obs[cluster_key].map(mapping).map({'False':False,'True':True})

    if override_info is not None:
        assert isinstance(override_info, dict) and 'col' in override_info and 'name' in override_info
        col,name = override_info['col'],override_info['name']
        assert isinstance(adata.obs[col],pd.Series) or isinstance(adata.obs[col],pd.Categorical)
        if name not in adata.obs[col].cat.categories:
            adata.obs[col] = adata.obs[col].cat.add_categories(name)
        unique_overrides = adata[adata.obs['is_doublet']].obs[col].cat.categories
        print(f"Replacing the values {unique_overrides} in column {col} with {name}")
        adata.obs.loc[adata.obs['is_doublet'],col] = name


    return mapping



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


def annotate_adata(adata,groups,other_group='?',cluster_col='leiden',class_col='class',output_pre='NR_',override=None,plot=True,output_f=""):
    df_l = []
    g_lst=[]
    var_names = []
    for group in groups:
        g_lst.append(group)
        info = groups[group]
        if 'z-scores' in info:
            assert len(info['marker_genes']) == len(info['z-scores'])
            var_names.extend(info['marker_genes'])
            cutoffs = np.array(adata[:,info['marker_genes']].X) > np.array(info['z-scores'])

            adata.obs[group]=np.any(cutoffs,axis=1)
            agg = adata.obs.groupby(cluster_col,observed=False)[group].mean().reset_index(name=group)
            agg.set_index(cluster_col,inplace=True)
            if plot:
                plt.hist(agg[group])
                plt.ylabel('Number of clusters')
                plt.xlabel('Prop of cells above cutoff zscore for at least 1 {} marker gene'.format(group))
                plt.show()

            potential_clusters = agg[agg[group] >info['pct_in']].index.values
            adata.obs[group] = pd.Categorical(adata.obs[cluster_col].isin(potential_clusters),categories=[False,True])

            if plot:
                sc.pl.umap(adata,color=[cluster_col]+info['marker_genes'] + [group],vmin=0,vmax=2)

        df_l.append(agg)
    if override is not None:
        for clust in override:
            ct = override[clust]
            if ct not in g_lst:
                g_lst.append(ct)
                adata.obs[ct] = False

            other_groups = list(set(g_lst) - {ct})
            adata.obs.loc[adata.obs[cluster_col]==clust,other_groups]=False
            adata.obs.loc[adata.obs[cluster_col]==clust, ct] = True
            adata.obs[ct] = pd.Categorical(adata.obs[ct],categories=[False,True])


    #return adata.obs[['leiden'] + g_lst]
    bool_t_ident = adata.obs[g_lst].values
    if np.any(np.sum(bool_t_ident,axis=1) >= 2):
        print("One cluster is assigned to two classes! Fix this. (error cells returned to df)")
        return adata.obs.iloc[np.nonzero(bool_t_ident.sum(axis=1)>=2)]
    #return adata.obs[['leiden'] + g_lst]

    assert not np.any(np.sum(bool_t_ident,axis=1) >= 2),"One cluster is assigned to two classes"
    
    idx_identities = ( np.sum(bool_t_ident,axis=1) !=0 ) * (bool_t_ident.argmax(axis=1)+1) #makes 0 the none type, 1 excitatory, 2 inhibitory...
    mapping = dict( zip(range(len(g_lst)+1), [other_group] + g_lst ) )
    print(mapping)
    adata.obs[class_col] = pd.Categorical( pd.Series(idx_identities,index=adata.obs.index).apply(lambda x: mapping[x]) )
    joint_col = f'{class_col}_{cluster_col}'
    adata.obs[joint_col] = adata.obs[class_col].astype(str) + "_" + adata.obs[cluster_col].astype(str)
    sc.pl.dotplot(adata,groupby=joint_col,var_names=var_names,vmin=0,vmax=2,use_raw=False)
    del adata.obs[joint_col]
    if output_f != "":
        adata.write_h5ad(output_f)

    

    
    return pd.concat(df_l,axis=1) 
    
def split_and_save(adata,col,data_dir,categories=None,detect_doublets=True,override_info=None,output_pre=""):
    if categories is None:
        categories=adata.obs[col].cat.categories
    for cat in categories:
        adata_c = adata[adata.obs[col]==cat].copy()
        adata_c.X = adata_c.raw.X
        del adata_c.uns,adata_c.obsm,adata_c.obsp,adata_c.varm
        pipeline(
            adata_c,batch=False,key_added='subleiden',output_f= os.path.join(data_dir,"{}{}.h5ad".format(output_pre,cat)),
            make_qc=False,detect_doublets=detect_doublets,override_info=override_info,inplace=True
        )
        sc.pl.umap(adata_c,color='subleiden')

def postprocess_and_merge(adata_merged,adata_lst,subclass_col,subtype_col,rep_col,subclass_order,cluster_key,detect_doublets,override_info={'col':'subclass','name':'Doublet'},output_f=""):
    adata_merged.obs[subtype_col] = pd.NA
    adata_merged.obs[subclass_col] = pd.NA
    
    for adata_f in adata_lst:
        assert isinstance(adata_f,str)
        adata = sc.read_h5ad(adata_f)
        if detect_doublets:
            run_scrublet(adata,rep_col=rep_col,cluster_key=cluster_key,recalc=False, score_cutoff=None,override_info=override_info)
            adata.write_h5ad(adata_f)
        class_subtypes = sue.create_subtypes(adata,groupby=subclass_col,key=cluster_key,subgroups=subclass_order,key_added=subtype_col)
        print(adata.obs[subtype_col].value_counts())
        adata_merged.obs.update(adata.obs)
        adata_merged.write_h5ad(output_f)

def merge_and_reprocess(adata_lst,joint_cluster_key_added,output_f,figure_dir,post_merge_func=None):
    to_merge=[]
    for adata_f in adata_lst:
        assert isinstance(adata_f,str)
        adata=sc.read_h5ad(adata_f)
        to_merge.append(adata)

    adata = ad.concat(to_merge,join='outer')
    if post_merge_func is not None:
        post_merge_func(adata)
   
    pipeline(adata,batch=False,min_genes=100,key_added=joint_cluster_key_added,output_f=output_f,figure_dir=figure_dir,inplace=True,make_qc=False)
    return adata

                    




    

