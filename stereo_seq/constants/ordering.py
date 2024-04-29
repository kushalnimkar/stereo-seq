
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

EXCLUDED_GROUPS = ['Doublet', '?']

def order_ctypes(ctypes, separator="_"):
    """
    Orders a list of ctypes based on a predefined class order and a secondary numerical value.

    Parameters:
        ctypes (list of str): List of ctype strings to be ordered.
        separator (str): The separator used to split the ctype into primary class and secondary value.

    Returns:
        list: A sorted list of ctypes based on primary class and secondary numerical value.
    """
    
    # Define a dictionary to hold secondary values extracted from ctypes
    # If the ctype does not contain the separator, assign a default value of -1
    secondary_keys = {
        ctype: int(ctype.split(separator)[1]) if separator in ctype else -1
        for ctype in ctypes
    }

    # Function to extract sort keys for each ctype
    def sort_key(ctype):
        primary_class, *_ = ctype.split(separator)
        primary_index = CLASS_ORDER.index(primary_class)
        secondary_value = secondary_keys[ctype]
        return (primary_index, secondary_value)

    # Return the list of ctypes sorted by the primary class order and secondary value
    return sorted(ctypes, key=sort_key)


def ordered_marker_genes(ordered_ctypes, separator="_"):
    """
    Generates a list of marker genes based on ordered ctypes and excluding specified groups.

    Parameters:
        ordered_ctypes (list of str): List of ordered ctypes taken from the function order_ctypes which takes unique values from the column 'subtypes'.
        separator (str): Separator used to split the ctype into primary class and secondary value.

    Returns:
        list: List of marker genes for the remaining ctypes after excluding specified groups.
    """

    # Extract unique primary classes from ctypes
    unique_classes = list(set([ctype.split(separator)[0] for ctype in ordered_ctypes]))

    # Filter out classes that are in the EXCLUDED_GROUPS
    included_classes = [c_class for c_class in unique_classes if c_class not in EXCLUDED_GROUPS]

    # Ensure that all remaining classes are in the CLASS_ORDER list
    assert all([c_class in CLASS_ORDER for c_class in included_classes]), f"One of {included_classes} is not in CLASS_ORDER"

    # Sort the included classes based on their index in CLASS_ORDER
    included_classes.sort(key=lambda x: CLASS_ORDER.index(x))

    # Initialize the list to hold marker genes
    marker_genes = []

    # Append marker genes from each class that is included
    for c_class in included_classes:
        if scm.MARKER_GENES[c_class]['marker_genes']:  # Ensure there are marker genes to append
            marker_genes.extend(scm.MARKER_GENES[c_class]['marker_genes'])

    return marker_genes



