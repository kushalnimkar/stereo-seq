MARKER_GENES={
    "Excitatory Neuron":{
        'marker_genes': ['Slc17a6'],
        'z-scores':[0.7],
        'pct_in': 0.2,
    },
    "Inhibitory Neuron":{
        'marker_genes':['Gad1','Gad2'],
        'z-scores':[1.0,1.0],
        'pct_in':0.4,
    },
    'Astrocyte':{
        'marker_genes':['Aqp4','Aldh1a1'],
        'z-scores':[1,2],
        'pct_in':0.4
    },
    'Microglia': {
        'marker_genes':['Cx3cr1'],#CCL3
        'z-scores':[2],
        'pct_in':0.4
    },
    'OPC':{
        'marker_genes':['Pdgfra'],
        'z-scores':[2],
        'pct_in':0.5,
    },
    'Oligodendrocyte':{
        'marker_genes':['Mbp'],
        'z-scores':[1.5],
        'pct_in':0.4
    },
    'Pericyte':{
        'marker_genes':['Cspg4','Pdgfrb'],
        'z-scores':[2,2],
        'pct_in':0.5,
    },
    'Meninges':{
        'marker_genes':['Foxc1'],
        'z-scores':[1.5],
        'pct_in':0.5,
    },
    'Epithelial':{
        'marker_genes':['Ranbp3l','Col4a5','Col4a6'],
        'z-scores':[2,2,2],
        'pct_in':0.4,
    },
    'T-cell':{
        'marker_genes':['Mrc1','F13a1'],
        'z-scores':[2,2],
        'pct_in':0.3
    },
    'Endothelial':{
        'marker_genes':['Slco1a4', 'Flt1'],
        'z-scores':[2,2],
        'pct_in':0.3
    }
}

def get_marker_genes():
    var_names = []
    [var_names.extend(MARKER_GENES[x]['marker_genes']) for x in MARKER_GENES]
    return var_names