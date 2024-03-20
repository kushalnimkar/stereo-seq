UNCERTAIN_MAPPING={
    'SC':{
        'P42':{
            'NR':{
                'clust_class_override':{'13':'Oligodendrocyte','29':'Pericyte','30':'Oligodendrocyte','31':'Oligodendrocyte'},
                'Non-neuronal_override':{
                    '0':'Oligodendrocyte','2':'Oligodendrocyte','6':'Oligodendrocyte','7':'Oligodendrocyte',
                    '3':'Doublet','8':'Doublet','11':'Doublet','12':'Doublet',
                    '13':'Meninges',
                    '14':'Pericyte',
                },
                # 'Non-neuronal_override':{
                #     '0':'Oligodendrocyte','11':'OPC','22':'Oligodendrocyte','26':'Doublet','28':'Oligodendrocyte','29':'Pericyte',
                #     '30':'Doublet','31':'Doublet','13':'Doublet'
                # },
            },
            'DR':{
                'clust_class_override':{'13':'Oligodendrocyte','15':'OPC','25':'Astrocyte','26':'T-cell','27':'Oligodendrocyte','28':'Pericyte'},
                'Non-neuronal_override':{
                    '0':'Oligodendrocyte','2':'Oligodendrocyte','3':'Oligodendrocyte','8':'Oligodendrocyte','10':'Oligodendrocyte',
                    '1':'Astrocyte',
                    '4':'OPC',
                    '10':'Meninges',
                    '11':'Doublet',
                    '13':'Doublet',
                    '5':'Doublet',
                    '6':'Doublet'
                },
            }
        } ,
    }
}

NR_NN_MAP = {
    '0':'Astrocyte',
    '1':'Oligodendrocyte',
    '2':'Oligodendrocyte',
    '3':'Oligodendrocyte',
    '4':'Doublet',
    '5':'Microglia',
    '6':'OPC',
    '7':'Oligodendrocyte',
    '8':'Doublet',
    '9':'Epithelial',
    '10':'Oligodendrocyte',
    '11':'Doublet',
    '12':'Doublet',
    '13':'Meninges',
    '14':'Pericyte',
    '15':'T-cell'
}

DR_NN_MAP={
    '0':'Astrocyte',
    '1':'Oligodendrocyte',
    '2':'Oligodendrocyte',
    '3':'Oligodendrocyte',
    '4':'OPC',
    '5':'Doublet',
    '6':'Doublet',
    '7':'Oligodendrocyte',
    '8':'Microglia',
    '9':'Oligodendrocyte',
    '10':'Epithelial',
    '11':'T-cell',
    '12':'Pericyte',#Doublet
    '13':'Meninges',
    '14':'Doublet',#Oligo/Astro doublet
}

DR_UNCERTAIN_MAPPING = {
    '6':'Astrocyte',
    '25':'Astrocyte',#Seems like neuron marker or doublet?
    '26':'Astrocyte',
    '27':'Microglia',#M2 microglia marker or OLIGO since it has mbp
    '28':'Pericyte' ,#Pdgfrb+, also seems to have some FTL1 (endo?)
}

NR_UNCERTAIN_MAPPING= {
    '4': 'Astrocyte',#Aldha1
    '13':'Oligodendrocyte',#Half excitatory, half inhibitory, but all have plp1 and MBP at higher levels
    '24':'Astrocyte',#Aldha1 (OR MENINGES)
    '25':'Microglia', #M2 microglia or macrophage, maybe inhibitory
    '26':'Astrocyte', #OR astrocyte+inhibitory doublet
    '28':'Pericyte', #Pdgfrb
    '29':'Oligodendrocyte',#Mbp and Aldha1+ positive Astrocyte or oligodendrocyte,
    '30':'Oligodendrocyte' #Could be microglia too
}

SUBClASS_TO_CLASS_MAPPING ={
    'Astrocyte':'Non-neuronal',
    'Excitatory Neuron':'Excitatory Neuron',
    'Inhibitory Neuron':'Inhibitory Neuron',
    'Microglia':'Non-neuronal',
    'OPC':'Non-neuronal',
    'Oligodendrocyte':'Non-neuronal',
    'Pericyte':'Non-neuronal',
    'Endothelial':'Non-neuronal',
    'Epithelial':'Non-neuronal',
    'Meninges':'Non-neuronal',
    'T-cell':'Non-neuronal'
}

CLASSES = sorted(list(set(SUBClASS_TO_CLASS_MAPPING.values())))
