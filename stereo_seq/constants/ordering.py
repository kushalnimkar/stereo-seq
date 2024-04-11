
import stereo_seq.constants.marker_genes as scm
NR_ORDERINGS ={
    'class':{
        'order':['Slc17a6','Gad1','Gad2','Aqp4','Aldh1a1','Cx3cr1','Pdgfra','Mbp','Pdgfrb','Cldn5','Foxc1','Ccdc146','Ranbp3l','Col4a5','Col4a6','Mrc1','F13a1']
    },
    'final_subclass':{
        'order':[
            *[f'Excitatory Neuron_{i}' for i in range(1,10)],
            *[f'Inhibitory Neuron_{i}' for i in range(1,11)],
            'Astrocyte',
            'Microglia',
            'OPC',
            *[f'Oligodendrocyte_{i}' for i in range (1,6)],
            'Pericyte',
            'Meninges',
            'Epithelial',
            'T-cell',
            *[f'Doublet_{i}' for i in range (1,5)]
        ],
        'marker_genes':None
    }
}

DR_ORDERINGS ={
    'final_subclass':{
        'order':[
            *[f'Excitatory Neuron_{i}' for i in range(1,10)],
            *[f'Inhibitory Neuron_{i}' for i in range(1,9)],
            'Astrocyte',
            'Microglia',
            'OPC',
            *[f'Oligodendrocyte_{i}' for i in range (1,6)],
            'Pericyte',
            'Meninges',
            'Epithelial',
            'T-cell',
            *[f'Doublet_{i}' for i in range (1,4)]
        ],
        'marker_genes':None      
    }
}

CLASS_ORDER = ['Excitatory Neuron', 'Inhibitory Neuron', 'Astrocyte','Microglia','OPC','Oligodendrocyte','Pericyte','Meninges','Epithelial','T-cell', 'Endothelial', 'Doublet', '?']

def order_ctypes(ctypes,separator="_"):
    secondary_keys = {x:int(x.split(separator)[1]) if separator in x else -1  for x in ctypes}
    return sorted(ctypes, key= lambda x: ( CLASS_ORDER.index(x.split(separator)[0]), secondary_keys[x]))


def ordered_marker_genes(ordered_ctypes,separator="_"):
    unique_classes = list(set([c_type.split(separator)[0] for c_type in ordered_ctypes]))
    assert all([c_class in CLASS_ORDER for c_class in unique_classes]), print(f"One of {unique_classes} is not in CLASS_ORDER onstant")
    unique_classes.sort(key= lambda x: CLASS_ORDER.index(x))
    marker_genes=[]
    [marker_genes.extend(scm.MARKER_GENES[c_class]['marker_genes']) for c_class in unique_classes]

    return marker_genes


