import numpy as np
from matplotlib import gridspec
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from sklearn.metrics.cluster import adjusted_rand_score
from collections import OrderedDict

PLT_CONTEXT={"figure.dpi":300,"font.size":12,'font.sans-serif':'Arial'}
PLT_SMALLER_CONTEXT={"figure.dpi":300,"font.size":9,'font.sans-serif':'Arial'}


def crosstab_confusion(adata,x_col,y_col,xlabel='',ylabel='',output_f=""):
    x_cats = adata.obs[x_col].cat.categories
    x_mapping = dict(zip(x_cats, range(len(x_cats))))

    y_cats= adata.obs[y_col].cat.categories
    y_mapping = dict(zip(y_cats, range(len(y_cats))))
    plot_mapping_new(
        test_labels=adata.obs[y_col].map(y_mapping).values.astype(int), #y-axis
        test_predlabels=adata.obs[x_col].map(x_mapping).values.astype(int),#x-AXIS
        test_dict=y_mapping,  # y-axis
        train_dict=x_mapping, # x-axis
        re_order=True,
        re_order_cols=None,
        re_index=False,
        xaxislabel=xlabel, yaxislabel=ylabel,
        save_as=output_f
)

def plotConfusionMatrixNew(
    ytrue,
    ypred,
    type,
    xaxislabel,
    yaxislabel,
    title,
    train_dict,
    test_dict=None,
    re_order=None,
    re_order_cols = None,
    re_index = None,
    re_order_rows = None,
    save_as=None,
    c='blue',
    ):
    # Note: Train classes = actual classes + unassigned
    # ytrue = independent clustering (test), ypred = classifier prediction (train)
    # Class labels must be nonnegative consecutive integers ranging from 0,1,...,n_classes-1


    numbertrainclasses = len(train_dict)
    numbertestclasses = len(test_dict)
    confusion_matrix = np.zeros((numbertestclasses,numbertrainclasses))
    for i,_ in enumerate(ytrue):
        confusion_matrix[ytrue[i],ypred[i]] +=1

    # Normalize confusion matrix
    for row,arr in enumerate(confusion_matrix):
        row_sum = np.sum(arr)
        if row_sum !=0:
            confusion_matrix[row] = confusion_matrix[row]/row_sum

    conf_df = pd.DataFrame(confusion_matrix)

    conf_df.index = list(test_dict.keys()) 
    conf_df.columns = list(train_dict.keys())

    # Reorder rows to try and make diagonal
    if re_index:
        raise NotImplementedError("under construction")
        inv_dict = {v:k for k,v in test_dict.items()}
        most_likely = np.argmax(confusion_matrix, axis=0)
        row_order = [inv_dict[i] for i in list(set(most_likely))]
        print(row_order)
        #row_order = list(dict.fromkeys(most_likely)) # Note: If one type is most likely for multiple rows, this will get rid of the duplicates
        unclear_assignment = set(test_dict.keys()) - set(row_order)
        row_order.extend(unclear_assignment)
        #row_order = [inv_dict[i] for i in row_order]
        print(row_order)
        conf_df = conf_df.reindex(row_order)
        print("Reindexing", row_order)

    if re_order and re_order_cols is None:
        most_likely_col_indices = np.argmax(conf_df.values, axis=1)
        most_likely_col_indices = list(dict.fromkeys(most_likely_col_indices))
        inv_dict = {v:k for k,v in train_dict.items()}
        col_order = [inv_dict[i] for i in most_likely_col_indices]
        #row_order = list(dict.fromkeys(most_likely)) # Note: If one type is most likely for multiple rows, this will get rid of the duplicates
        unclear_assignment = set(train_dict.keys()) - set(col_order)
        col_order.extend(unclear_assignment)
        try:
            col_order.remove('Unassigned')
            col_order.append('Unassigned')
        except:
            pass
        conf_df = conf_df[col_order]
        print("Reiordering columns to try and make diagonal")
    elif re_order:        
        conf_df=conf_df[re_order_cols]
    
    diagcm = conf_df.to_numpy()
    xticksactual = conf_df.columns

    dot_max = np.max(diagcm.flatten())
    dot_min = 0
    if dot_min != 0 or dot_max != 1:
        frac = np.clip(diagcm, dot_min, dot_max)
        old_range = dot_max - dot_min
        frac = (frac - dot_min) / old_range
    else:
        frac = diagcm
    xvalues = []
    yvalues = []
    sizes = []
    for i in range(diagcm.shape[0]):
        for j in range(diagcm.shape[1]):
            xvalues.append(j)
            yvalues.append(i)
            sizes.append((frac[i,j]*35)**1.5)
    size_legend_width = 0.5
    height = diagcm.shape[0] * 0.3 + 1
    height = max([1.5, height])
    heatmap_width = diagcm.shape[1] * 0.35
    width = (
        heatmap_width
        + size_legend_width
        )
    fig = plt.figure(figsize=(width, height))
    axs = gridspec.GridSpec(
        nrows=2,
        ncols=2,
        wspace=0.02,
        hspace=0.04,
        width_ratios=[
                    heatmap_width,
                    size_legend_width
                    ],
        height_ratios = [0.5, 10]
        )
    dot_ax = fig.add_subplot(axs[1, 0])
    dot_ax.scatter(xvalues,yvalues, s = sizes, c = c, norm=None, edgecolor='none')
    y_ticks = range(diagcm.shape[0])
    dot_ax.set_yticks(y_ticks)
    if type == 'validation':
        dot_ax.set_yticklabels(list(train_dict.keys()))
    elif type == 'mapping':
        # dot_ax.set_yticklabels(list(test_dict.keys()))
        dot_ax.set_yticklabels(list(conf_df.index))
    
    #xticksactual = conf_df.columns.map({v:k for k,v in train_dict.items()})
    dot_ax.set_xticklabels(xticksactual,rotation=90)
    x_ticks = range(diagcm.shape[1])
    dot_ax.set_xticks(x_ticks)
    #dot_ax.set_xticklabels(xticksactual, rotation=90)
    dot_ax.tick_params(axis='both', labelsize='small')
    dot_ax.grid(True, linewidth = 0.2)
    dot_ax.set_axisbelow(True)
    dot_ax.set_xlim(-0.5, diagcm.shape[1] + 0.5)
    ymin, ymax = dot_ax.get_ylim()
    dot_ax.set_ylim(ymax + 0.5, ymin - 0.5)
    dot_ax.set_xlim(-1, diagcm.shape[1])
    dot_ax.set_xlabel(xaxislabel)
    dot_ax.set_ylabel(yaxislabel)
    dot_ax.set_title(title)
    size_legend_height = min(1.75, height)
    wspace = 10.5 / width
    axs3 = gridspec.GridSpecFromSubplotSpec(
        2,
        1,
        subplot_spec=axs[1, 1],
        wspace=wspace,
        height_ratios=[
                    size_legend_height / height,
                    (height - size_legend_height) / height
                    ]
        )
    diff = dot_max - dot_min
    if 0.3 < diff <= 0.6:
        step = 0.1
    elif diff <= 0.3:
        step = 0.05
    else:
        step = 0.2
    fracs_legends = np.arange(dot_max, dot_min, step * -1)[::-1]
    if dot_min != 0 or dot_max != 1:
        fracs_values = (fracs_legends - dot_min) / old_range
    else:
        fracs_values = fracs_legends
    size = (fracs_values * 35) ** 1.5
    size_legend = fig.add_subplot(axs3[0])
    size_legend.scatter(np.repeat(0, len(size)), range(len(size)), s=size, c = c)
    size_legend.set_yticks(range(len(size)))
    labels = ["{:.0%}".format(x) for x in fracs_legends]
    if dot_max < 1:
        labels[-1] = ">" + labels[-1]
    size_legend.set_yticklabels(labels)
    size_legend.set_yticklabels(["{:.0%}".format(x) for x in fracs_legends])
    size_legend.tick_params(axis='y', left=False, labelleft=False, labelright=True)
    size_legend.tick_params(axis='x', bottom=False, labelbottom=False)
    size_legend.spines['right'].set_visible(False)
    size_legend.spines['top'].set_visible(False)
    size_legend.spines['left'].set_visible(False)
    size_legend.spines['bottom'].set_visible(False)
    size_legend.grid(False)
    ymin, ymax = size_legend.get_ylim()
    size_legend.set_ylim(ymin, ymax + 0.5)
    if save_as is not None:
        fig.savefig(save_as, bbox_inches = 'tight')
    plt.show()
    return diagcm, None,axs

def plot_mapping_new(test_labels, test_predlabels, test_dict, train_dict, xaxislabel, yaxislabel,re_order=None,re_order_cols = None,re_index = None, re_order_rows = None, save_as=None,c='blue'):
    
    ARI = adjusted_rand_score(labels_true = test_labels, 
                              labels_pred = test_predlabels)
    accuracy = np.sum(test_labels==test_predlabels)/len(test_labels)*100
    
           
    mappingconfmat, mappingxticks, mappingplot = plotConfusionMatrixNew(
    ytrue = test_labels,
    ypred = test_predlabels,
    test_dict=test_dict,
    train_dict=train_dict,
    type = 'mapping',
    save_as = save_as,
    title = 'ARI = {:.3f}, Accuracy = {:.1f}%'.format(ARI,accuracy),
    xaxislabel =xaxislabel,
    yaxislabel = yaxislabel,
        re_order=re_order,
    re_order_cols = re_order_cols,
        re_index = re_index,
    re_order_rows = re_order_rows,
    c=c,
    )

def construct_dotplot(adata,marker_genes,cluster_key,exclude=None,include=None):
    for crit in [include,exclude]:
        if crit is not None:
            assert all([x in adata.obs[cluster_key].cat.categories for x in crit])
    
    
    if exclude is not None:
        cats = [c for c in adata.obs[cluster_key].cat.categories if c not in exclude]
    elif include is not None:
        cats=include
    else:
        cats = adata.obs[cluster_key].cat.categories

    var_group_labels = [marker_genes[c] for c in cats]
    var_names = sum(var_group_labels,[])
    #print(var_names)
    var_names = list(dict.fromkeys(var_names))
    print(var_names)
    adata=adata[adata.obs[cluster_key].isin(cats)]
    sc.pl.dotplot(adata,groupby=cluster_key,var_names=sum(var_group_labels,[]))
