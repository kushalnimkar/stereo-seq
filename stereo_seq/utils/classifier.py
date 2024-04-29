import os
import numpy as np
import pandas as pd
import xgboost as xgb
import pickle as pickle
import time
import stereo_seq.utils.plotting as pl
import stereo_seq.utils.classifier as cl
import stereo_seq.utils.preprocessing as pre
import stereo_seq.models.model_info as smi

import importlib

def train_model(adata,obs_id, train_dict, model_dir,model_name, cutoff, common_top_genes=None, eta=0.2,max_cells_per_ident=2000,train_frac=0.7,min_cells_per_ident=100,run_validation=True,run_full=False,unassigned_key="Unassigned",colsample_bytree=1,nround=200,max_depth=4):
    if common_top_genes is None:
        common_top_genes= list(adata.var.index)

    y, y_probs, valid_model,model=trainclassifier(
        train_anndata=adata, 
        common_top_genes=common_top_genes,                                                 
        obs_id=obs_id, # Clusters
        train_dict=train_dict,  #Cluster: value
        eta=eta,
        max_cells_per_ident=max_cells_per_ident, train_frac=train_frac,
        min_cells_per_ident=min_cells_per_ident,
        run_validation=run_validation,
        run_full=run_full,
        colsample_bytree=colsample_bytree,
        nround=nround,
        max_depth=max_depth
    )

    if model is not None:
        model.save_model(os.path.join(model_dir,model_name))
    if valid_model is not None:
        valid_model.save_model(os.path.join(model_dir,f"valid_{model_name}"))

    y_pred = threshold(y_probs,unassigned_indicator=train_dict[unassigned_key],cutoff=cutoff)                                            
    np.save(os.path.join(model_dir,f"{model_name}_valid_y"),y)
    np.save(os.path.join(model_dir,f"{model_name}_valid_y_pred"),y_pred)
    if model is not None:
        model.save_model(os.path.join(model_dir,model_name))
    if valid_model is not None:
        valid_model.save_model(os.path.join(model_dir,f"valid_{model_name}"))
    return y,y_pred,valid_model,model


def gen_dict(adata,col,prefix,unassigned_indicator):
    classes = adata.obs[col].unique().categories

    if all([i.isnumeric() for i in classes]):
        classes = sorted(classes.astype(int))

    if prefix is not None:
        classes = [f'{prefix}{i}' for i in classes]

    n_classes = len(classes)
    t_dict= dict(zip(classes, range(n_classes)))

    if unassigned_indicator is not None:
        t_dict[unassigned_indicator]=n_classes

    return t_dict
    
    
    


#This helper method uses xgboost to train classifiers.
def trainclassifier(train_anndata, common_top_genes, obs_id, train_dict, eta,
                    max_cells_per_ident, train_frac, min_cells_per_ident, run_validation,run_full,colsample_bytree=1,nround=200,max_depth=4):
    
    start_time = time.time()
    unique_classes = np.unique(train_anndata.obs[obs_id].values)
    numbertrainclasses = len(unique_classes)

    xgb_params_train = {
            'objective':'multi:softprob',
            'eval_metric':'mlogloss',
            'num_class':numbertrainclasses,
            'eta':eta,
            'max_depth':max_depth,
            'subsample': 0.6,
            'colsample_bytree':colsample_bytree,
            'n_jobs': 6}
    #Train XGBoost on 70% of training data and validate on the remaining data


    training_set_train_70 = []
    validation_set_train_70 = []
    training_label_train_70 = []
    validation_label_train_70 = []
    bst_model_train_70 =  None
    bst_model_full_train = None
    validation_pred_train_70 = None
    if run_validation:
        print("Running train/validation split")
        #loop thru classes to split for training and validation
        for i in unique_classes:

            #how many cells in a class
            cells_in_clust = train_anndata[train_anndata.obs[obs_id]==i,:].obs_names #cell names
            n = min(max_cells_per_ident,round(len(cells_in_clust)*train_frac))

            #sample 70% for training and rest for validation
            train_temp = np.random.choice(cells_in_clust,n,replace = False)
            validation_temp = np.setdiff1d(cells_in_clust, train_temp)

            #upsample small clusters
            if len(train_temp) < min_cells_per_ident:
                train_temp_bootstrap = np.random.choice(train_temp, size = min_cells_per_ident - int(len(train_temp)))
                train_temp = np.hstack([train_temp_bootstrap, train_temp])

            #store training and validation **names** of cells in vectors, which update for every class
            training_set_train_70 = np.hstack([training_set_train_70,train_temp])
            validation_set_train_70 = np.hstack([validation_set_train_70,validation_temp])

            #store training and validation **labels** of cells in vectors, which update for every class
            training_label_train_70 = np.hstack([training_label_train_70,np.repeat(train_dict[i],len(train_temp))])
            validation_label_train_70 = np.hstack([validation_label_train_70,np.repeat(train_dict[i],len(validation_temp))])

            #need train_dict b/c XGboost needs number as class labels, not words
            #this is only deconvulted later in plotting function

        #put data in XGB format
        X_train = train_anndata[training_set_train_70,common_top_genes].X
        train_matrix_train_70 = xgb.DMatrix(data = X_train, label = training_label_train_70, 
                                            feature_names = common_top_genes)

        X_valid = train_anndata[validation_set_train_70,common_top_genes].X
        validation_matrix_train_70 = xgb.DMatrix(data = X_valid, label = validation_label_train_70, 
                                                 feature_names = common_top_genes)

        del training_set_train_70, validation_set_train_70, training_label_train_70

        #Train on 70%
        bst_model_train_70 = xgb.train(
            params = xgb_params_train,
            dtrain = train_matrix_train_70,
            num_boost_round = nround)
        
        validation_pred_train_70 = bst_model_train_70.predict(data = validation_matrix_train_70)
    if run_full:
        print("Training on 100%")
        #Train on 100%
        #Train XGBoost on the full training data
        training_set_train_full = []
        training_label_train_full = []

        for i in unique_classes:
            train_temp = train_anndata.obs.index[train_anndata.obs[obs_id].values == i]
            if len(train_temp) < 100:
                train_temp_bootstrap = np.random.choice(train_temp, size = 100 - int(len(train_temp)))
                train_temp = np.hstack([train_temp_bootstrap, train_temp])
            
            #indices of cells in class
            training_set_train_full = np.hstack([training_set_train_full,train_temp])
            
            #labels of cells in class: [label*N_class] stacked onto previous classes
            training_label_train_full = np.hstack([training_label_train_full,np.repeat(train_dict[i],len(train_temp))])


        X_train_full = train_anndata[training_set_train_full,common_top_genes].X
        full_training_data = xgb.DMatrix(data = X_train_full, label = training_label_train_full, 
                                        feature_names = common_top_genes)

        del training_set_train_full, training_label_train_full

        bst_model_full_train = xgb.train(
            params = xgb_params_train,
            dtrain = full_training_data,
            num_boost_round = nround
            )

    
    
    print('trainclassifier() complete after', np.round(time.time() - start_time), 'seconds')
    
    #real labels of validation set, predicted labels, classifier.
    #recall these are all integers that are deconvulted later in plotting using the dicts
    return validation_label_train_70, validation_pred_train_70, bst_model_train_70, bst_model_full_train

"""Vignette"""
def run_training_and_validation(adata,obs_id, model_dir,model_name,figure_dir,output_f, train_dict_name, run_validation=True,run_full=False,xlabel='Predicted',ylabel='True',
                cutoff=0.8,eta=0.2,max_cells_per_ident=2000,train_frac=0.7,min_cells_per_ident=100, colsample_bytree=1,nround=200, max_depth=4, unassigned_indicator='Unassigned'):
    
    train_dict = gen_dict(adata,obs_id,prefix='',unassigned_indicator=unassigned_indicator)
    print(train_dict)
    with open(os.path.join(model_dir,train_dict_name),'wb') as f:
        pickle.dump(train_dict,f)
    y,y_pred,valid_model,model = train_model(
        adata,
        obs_id, 
        train_dict, 
        model_dir,
        model_name, 
        cutoff, 
        common_top_genes=None, 
        eta=eta,
        max_cells_per_ident = max_cells_per_ident,
        train_frac=train_frac,
        min_cells_per_ident = min_cells_per_ident,
        run_validation = run_validation,
        run_full = run_full,
        unassigned_key=unassigned_indicator,
        colsample_bytree=colsample_bytree,
        nround=nround,
        max_depth=max_depth
    )
    #test_dict = gen_dict(adata,obs_id,prefix=None,unassigned_indicator=None)
    #print(train_dict,test_dict)
    test_dict = {k:v for k,v in train_dict.items() if k != unassigned_indicator}
    importlib.reload(pl)
    pl.plot_mapping_new(test_labels=y.astype(int), #y-axis
                            test_predlabels=y_pred,#x-AXIS
                            test_dict=test_dict,  # y-axis
                            train_dict=train_dict, # x-axis
                            re_order=True,
                            re_order_cols=None,
                            re_index=False,
                            xaxislabel=xlabel, yaxislabel=ylabel,
                            save_as=os.path.join(figure_dir,output_f)
    )
    return valid_model,model,train_dict,test_dict

def threshold(y_probs,unassigned_indicator,cutoff):
    maxes = np.max(y_probs,axis=1)
    max_indices = np.argmax(y_probs,axis=1)
    y_pred = np.where(maxes < cutoff, unassigned_indicator, max_indices)
    return y_pred


def apply_model(adata,model_f,train_dict,cutoff, figure_dir,output_f, unassigned_indicator='Unassigned', prefix='MC', test_key='leiden',key_added='retina_class',prefix_key=None,xlabel='Predicted',ylabel='Clustering'): #Test key is usually clustering
    if isinstance(train_dict,str):
        with open(train_dict,'rb') as f:
            train_dict = pickle.load(f)
    if isinstance(model_f, str):
        model=xgb.Booster()
        model.load_model(model_f)
    else:
        model=model_f # Otherwise directly provide model

    assert isinstance(train_dict,dict)
    test_dict = gen_dict(adata,col=test_key,prefix=prefix,unassigned_indicator=None)
    """
    if prefix_key is None:
        prefix_key = f"{prefix}_{test_key}"
    """


    y_pred,y_probs = cutoff_predict(X=adata.X, model=model, unassigned_indicator=train_dict[unassigned_indicator], cutoff=cutoff)
    test_labels = adata.obs[test_key]
    if adata.obs[test_key][0].isnumeric():
        test_labels=adata.obs[test_key].astype(int)
    else:
        test_labels = adata.obs[test_key].map(test_dict)

    adata.obs[key_added] = pd.Series(y_pred).map({v:k for k,v in train_dict.items()}).values
    print(test_dict,train_dict)
    pl.plot_mapping_new(test_labels=test_labels, #y-axis
                            test_predlabels=y_pred,#x-AXIS
                            test_dict=test_dict,  # y-axis
                            train_dict=train_dict, # x-axis
                            re_order=True,
                            re_order_cols=None,
                            re_index=False,
                            xaxislabel=xlabel, yaxislabel=ylabel,
                            save_as=os.path.join(figure_dir,output_f)
                )
    return test_dict

def cutoff_predict(X, model, unassigned_indicator, cutoff,gene_names=None):
    if gene_names is None:
        y_probs = model.predict(xgb.DMatrix(data = X))
    else:
        y_probs = model.predict(xgb.DMatrix(data = X,feature_names=gene_names))
    y_pred = threshold(y_probs,unassigned_indicator,cutoff)
    return y_pred,y_probs

def threshold(y_probs,unassigned_indicator,cutoff):
    maxes = np.max(y_probs,axis=1)
    max_indices = np.argmax(y_probs,axis=1)
    y_pred = np.where(maxes < cutoff, unassigned_indicator, max_indices)
    return y_pred



def train_and_test(adata_train,adata_test,train_col,train_subset, model_dir,model_name,figure_dir,train_dict_name,test_col=None,key_added=None,validation_output_f="validation.png",test_output_f="test_compare.png",genes=None,cutoff=0.5,run_full=False,preprocess_train=False,preprocess_test=False,train_model=True,evaluate_test=True):
    for d in [model_dir,figure_dir]:
        if not os.path.exists(d):
            os.mkdir(d)

    if genes is None:
        print("Note: Ensuring two datasets have same gene set. Make sure to rename any that have different names beforehand...")
        genes = sorted(list(set(adata_train.var.index).intersection(set(adata_test.var.index))))
            
    adata_train._inplace_subset_var(genes)
    adata_test._inplace_subset_var(genes)

    
    if train_subset is not None:
        adata_train._inplace_subset_obs(adata_train.obs[adata_train.obs[train_col].isin(train_subset)].index)
    

    if preprocess_train:
        print("Preprocessing training")
        pre.preprocess(adata=adata_train,min_genes=None,batch=False,max_counts=None,inplace=True,find_hvg=False,do_pca=False)
    if preprocess_test:
        print("Preprocessing test")
        pre.preprocess(adata=adata_test,min_genes=None,batch=False,max_counts=None,inplace=True,find_hvg=False,do_pca=False)
    
    adata_train.obs[train_col] = pd.Categorical(adata_train.obs[train_col])
    info={}
    if train_model:
        valid_model,model,train_dict,test_dict = cl.run_training_and_validation(
            adata=adata_train,
            obs_id=train_col,
            model_dir=model_dir, model_name=model_name,figure_dir=figure_dir,output_f =f"{cutoff}_validation.png",train_dict_name=train_dict_name,
            colsample_bytree=1,nround=200,max_cells_per_ident=1000,min_cells_per_ident=200,max_depth=4,
            run_full=run_full
        )
        info['valid_model'] = valid_model
        info['model'] = model
        info['train_dict'] = train_dict
        info['test_dict'] = test_dict
    
    if evaluate_test:
        if run_full:
            xgb_model_f=os.path.join(model_dir, f"{model_name}")
        else:
            xgb_model_f=os.path.join(model_dir, f"valid_{model_name}")
        with open(os.path.join(model_dir,train_dict_name),'rb') as f:
            train_dict= pickle.load(f)
            info['train_dict']=train_dict
        cl.apply_model(adata_test,xgb_model_f,train_dict,cutoff=cutoff, figure_dir=figure_dir,output_f=test_output_f, unassigned_indicator='Unassigned', prefix=None, test_key=test_col,key_added=key_added,prefix_key=None,xlabel='Predicted',ylabel='Previous Clustering') #Test key is usually clustering
    
    return info

def run_classification_pipeline(region,timepoint,model_name,test_data_dir,test_col,key_added,preprocess_train=True,preprocess_test=False,train_model=False,run_full=False,evaluate_test=True,resave=True):
    info = smi.get_model_info(region,timepoint,model_name)
    dataset_info = info['dataset']
    model_info = info['model']
    test_group=smi.MODEL_GROUP_MAP[model_name]#Consider replacing test_group and test_data_dir with test_f directly


    adata_train = dataset_info['loader'](dataset_info['path'])
    test_f=os.path.join(os.path.join(test_data_dir,f"{test_group}.h5ad"))
    adata_test = sc.read_h5ad(test_f)#Already preprocessed
    suc.train_and_test(
        adata_train,
        adata_test,
        **model_info,
        test_output_f="test_set.png",
        test_col=test_col,
        key_added=key_added,
        genes=None,
        cutoff=0.5,
        preprocess_train=preprocess_train,
        preprocess_test=preprocess_test,
        train_model=train_model,
        run_full=run_full,
        evaluate_test=evaluate_test,
    )

    if resave:
        adata_test.write_h5ad(test_f)
    return adata_test
    



    







