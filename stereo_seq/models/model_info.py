import os
import stereo_seq.constants.paths as scp
import stereo_seq.utils.load as sul

import scanpy as sc



_MODEL_INFO = {
    'SC':{
        'P42':{
            'excitatory_Xie_SC':{
                'dataset':{
                    'path':os.path.join(scp.REF_DATASETS,'Xie_SC.h5ad'),
                    'loader':sul.load_XieSC
                },
                'model':{
                    "model_name":"excitatory_Xie_SC",
                    'train_dict_name':"train_dict.pickle",
                    'validation_output_f':"validation.png",
                    'train_col':'celltype',
                    'train_subset':[f'Ex-{i}' for i in range(1,10)]
                }
            },
            'inhibitory_Xie_SC':{
                'dataset':{
                    'path':os.path.join(scp.REF_DATASETS,'Xie_SC.h5ad'),
                    'loader':sul.load_XieSC
                },
                'model':{
                    "model_name":"inhibitory_Xie_SC",
                    'train_dict_name':"train_dict.pickle",
                    'validation_output_f':"validation.png",
                    'train_col':'celltype',
                    'train_subset':[f'In-{i}' for i in range(1,11)]
                },
            },
            'excitatory_NR':{
                'dataset':{
                    'path':os.path.join(scp.make_path(base='data',region='SC',timepoint='P42',group='NR'),"Excitatory Neuron.h5ad"),
                    'loader':sc.read_h5ad
                },
                'model':{
                    "model_name":"excitatory_NR",
                    'train_dict_name':"train_dict.pickle",
                    'validation_output_f':"validation.png",
                    'train_col':'subleiden',
                    'train_subset':None
                }
            },
            'inhibitory_NR':{
                'dataset':{
                    'path':os.path.join(scp.make_path(base='data',region='SC',timepoint='P42',group='NR'),"Inhibitory Neuron.h5ad"),
                    'loader':sc.read_h5ad
                },
                'model':{
                    "model_name":"inhibitory_NR",
                    'train_dict_name':"train_dict.pickle",
                    'validation_output_f':"validation.png",
                    'train_col':'subleiden',
                    'train_subset':None
                }
            },
            'non-neuronal_NR':{
                'dataset':{
                    'path':os.path.join(scp.make_path(base='data',region='SC',timepoint='P42',group='NR'),"Non-neuronal.h5ad"),
                    'loader':sc.read_h5ad
                },
                'model':{
                    "model_name":"non-neuronal_NR",
                    'train_dict_name':"train_dict.pickle",
                    'validation_output_f':"validation.png",
                    'train_col':'subleiden',
                    'train_subset':None
                }
            }
        }
    }
}

def get_model_info(region,timepoint,model_name):
    info = _MODEL_INFO[region][timepoint][model_name]
    model_info = info['model']

    for key in ['model_dir','figure_dir']:
        if key not in model_info:
            model_info[key] = os.path.join(scp.MODEL_DIR,region,timepoint,model_name)
    
    if 'model_name' not in model_info:
        model_info['model_name'] = model_name
        
    
    for key in ['validation_output_f']:
        if key not in model_info:
            model_info[key] = "validation.png"
    
    if not os.path.exists(model_info['model_dir']):
        print(f"Making path for model since it does not exist: {model_info['model_dir']}")
        os.makedirs(model_info['model_dir'])
    
    if 'loader' not in info['dataset']:
        info['dataset']['loader'] = sc.read_h5ad

    return info


MODEL_GROUP_MAP={
    'excitatory_Xie_SC': 'Excitatory Neuron',
    'inhibitory_Xie_SC': 'Inhibitory Neuron',
    'excitatory_NR': 'Excitatory Neuron',
    'inhibitory_NR': 'Inhibitory Neuron',
    'non-neuronal_NR': 'Non-neuronal'
}
