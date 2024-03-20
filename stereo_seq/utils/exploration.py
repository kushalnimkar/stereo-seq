import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import scrublet as scr

import os
import re

import stereo_seq.constants.ordering as sco
from stereo_seq.utils.plotting import PLT_CONTEXT,PLT_SMALLER_CONTEXT

def gen_candidate_markers(adata,groupby,cats,n_genes=10,diff=True):
    if 'log1p' not in adata.layers:
        adata.layers['log1p'] = sc.pp.log1p(adata.layers['raw'])
    sc.tl.rank_genes_groups(adata,groupby=groupby,method='wilcoxon',layer='log1p',use_raw=False,pts=True)
    type_marker_map = dict()

    for n_type in cats:

        df = adata.uns['rank_genes_groups']['pts']

        if diff:
            pts_in =df[n_type]
            other_columns = list(set(df.columns) - {n_type})
            pts_out = np.max(df[other_columns].values,axis=1)
            pts_diff = pts_in-pts_out
            var_names=pts_diff.sort_values(ascending=False).index[0:n_genes]
        else:
            var_names = var_names=adata.uns['rank_genes_groups']['names'][n_type]

        with plt.rc_context({"font.sans-serif":"Arial"}):
            sc.pl.dotplot(adata,var_names=var_names,groupby=groupby,show=False,vmin=0,vmax=1)
            #plt.savefig("")
        type_marker_map[n_type] = list(var_names)

    return type_marker_map

def create_subtypes(adata,groupby,subgroups,key='leiden',key_added='class_subtypes',output_f=""):
    class_subgroups = adata.obs.groupby(groupby,observed=False)[key].unique()
    class_subtypes = []
    all_clusters=[]
    if subgroups is None:
        subgroups = class_subgroups.index.values

    for ct in subgroups:
        if ct not in class_subgroups:
            print(f"Warning, {ct} is not present as an actually annotated class; Skipping")
            continue
        leidens = sorted(class_subgroups[ct],key=int)
        all_clusters.extend(leidens)
        if len(leidens) == 1:
            post_added = [f"{ct}" for i in range(1,1+len(leidens))]
        else:
            post_added = [f"{ct}_{i}" for i in range(1,1+len(leidens))]
        class_subtypes.extend(post_added)

    missing_clusters = set(adata.obs[key].unique()) - set(all_clusters)
    if len(missing_clusters) != 0:
        missing_classes = set(adata.obs[groupby].unique()) - set(subgroups)
        raise AssertionError(f"Missing classes: {missing_classes}, missing_clusters = {missing_clusters}")
    mapping = dict(zip(all_clusters,class_subtypes))
    adata.obs[key_added] = adata.obs[key].map(mapping)
    return class_subtypes

def summary_stats(adata,subtype_col,color_col,var_names,figure_dir,rep_col='sample',categories_order=None,output_pre=""):
    if categories_order is None:
        categories_order = sco.order_ctypes(adata.obs[subtype_col].unique())
    if var_names is None:
        var_names= sco.ordered_marker_genes(categories_order)

        
    unique_subtypes= adata.obs[subtype_col].unique()
    for cat in categories_order:
        if cat not in unique_subtypes:
            print(f"Warning: Category in categories order: {cat} is not in unique col for {subtype_col}")
            categories_order.remove(cat)

    with plt.rc_context(PLT_CONTEXT):
        sc.pl.dotplot(adata,groupby=subtype_col,categories_order=categories_order,var_names=var_names,vmin=0,vmax=2,use_raw=False,show=False)
        plt.savefig(os.path.join(figure_dir,f"{output_pre}class_subgroups_dotplot.png"),bbox_inches='tight')
        plt.show()
    
        sc.pl.umap(adata,color=[rep_col,color_col],show=False)
        plt.savefig(os.path.join(figure_dir,f"{output_pre}umap.png"),bbox_inches='tight')
        plt.show()
    
    with plt.rc_context(PLT_CONTEXT):
        ax=sns.countplot(adata.obs,x=subtype_col,hue=color_col,stat='percent',order=adata.obs[subtype_col].value_counts().index)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
        ax.set_ylabel(f'% of cells')
        ax.set_title(f"Total cells: {len(adata.obs)} ({len(adata.obs[rep_col].unique())} replicates)")
        plt.savefig(os.path.join(figure_dir,f"{output_pre}proportions.png"),bbox_inches='tight')
        plt.show()
    
    with plt.rc_context(PLT_CONTEXT):
        ax=sns.countplot(adata.obs,x=color_col,hue=color_col,stat='percent',order=adata.obs[color_col].value_counts().index)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
        ax.set_ylabel(f'% of cells')
        ax.set_title(f"Total cells: {len(adata.obs)} ({len(adata.obs[rep_col].unique())} replicates)")
        plt.savefig(os.path.join(figure_dir,f"{output_pre}_class_only_proportions.png"),bbox_inches='tight')
        plt.show()
    
def summary_integration(adata,cluster_key,sample_key,figure_dir,pre="",output_data_f=""):
    with plt.rc_context(PLT_CONTEXT):
        sc.pl.umap(adata,color=[sample_key,cluster_key],show=False)
        plt.savefig(os.path.join(figure_dir,f"{pre}condition_leiden_umap.png"))
        plt.show()

        sns.countplot(adata.obs,x=cluster_key,hue=sample_key,stat='proportion')
        plt.xlabel('Joint NR and DR clusters')
        plt.ylabel('Proportion of cells')
        plt.savefig(os.path.join(figure_dir,f"{pre}cluster_proportions.png"))
        plt.show()

    adata.write_h5ad(output_data_f)
        

def gen_assignment(adata,col_0,col_1,assignment_strat=None,cutoff=0.5,unmapped_key='Unmapped'):#For every type in col_0, assign to col_1 (col_1 is usually classifier)
    cross = pd.crosstab(adata.obs[col_0],adata.obs[col_1])
    props = cross.values/np.sum(cross.values,axis=1)[:,None]
    unmapped_indicator = -1

    #Implement alternative assignment strategies

    if assignment_strat is None or assignment_strat=="max":
        assert cutoff <=1 and cutoff >=0
        indices = np.where(np.max(props,axis=1)>cutoff,np.argmax(props,axis=1),unmapped_indicator)#-1
    elif assignment_strat=="entropy":
        raise NotImplementedError(f"Assignment strategy: '{assignment_strat}' not yet implemented")
    
    mapped_type = np.where(indices != unmapped_indicator, np.array(cross.columns)[indices],'Unmapped')
    mapping = dict(zip(cross.index.values,mapped_type))



    return mapping

def get_stratified_diff_genes(adata,split_key,type_key,key_added,subcats=None,plot=False,pre="NN",output_dir=None):
    info = dict()
    info['type_key'] = type_key
    info['split_key'] = split_key
    info['diff_genes'] = {x:dict() for x in adata.obs[split_key].cat.categories}

    if subcats is None:
        subcats = adata.obs[type_key].cat.categories

    for cat in subcats:
        adata_t = adata[adata.obs[type_key]==cat]
        sc.tl.rank_genes_groups(adata_t,groupby=split_key,method='wilcoxon',layer='log1p',use_raw=False,pts=True)
        sc.tl.filter_rank_genes_groups(adata_t,key_added=key_added)
        if plot:
            df = sc.get.rank_genes_groups_df(adata_t,group=adata.obs[split_key].cat.categories,key='rank_genes_groups_filtered')
            if len(df) !=0:
                sc.pl.rank_genes_groups_dotplot(adata_t, key=key_added,vmin=0,vmax=1,use_raw=False,title=cat,show=False)
                if output_dir is not None:
                    plt.savefig(os.path.join(output_dir,f"{pre}{cat}.png"),bbox_inches='tight')
                    plt.show()
            else:
                print(f"Skipping {cat} since no diff expressed")
        for cat_group in adata_t.obs[split_key].cat.categories:
            df = sc.get.rank_genes_groups_df(adata_t,group=[cat_group],key='rank_genes_groups_filtered')

            if len(df) == 0:
                info['diff_genes'][cat_group][cat] = []
            else:
                info['diff_genes'][cat_group][cat] = list(df['names'].values)
    
    hits = _get_top_hits(info)
    return info,hits

def _get_top_hits(info):
    diff_genes = info['diff_genes']
    hits = dict()

    for cat in diff_genes.keys():
        hits[cat] = dict()
        for c in diff_genes[cat]:
            genes = diff_genes[cat][c]
            for g in genes:
                if g not in hits[cat]:
                    hits[cat][g] = 1
                else:
                    hits[cat][g] += 1
    return hits





