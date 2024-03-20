import os
BASE_DIR=".."
DATA_DIR,FIGURE_DIR,MODEL_DIR = os.path.join(BASE_DIR,"data"),os.path.join(BASE_DIR,"figures"),os.path.join(BASE_DIR,"models")
REF_DATASETS = os.path.join(DATA_DIR,"atlas")

BASES=['data','figures']
REGIONS=['SC']
TIMEPOINTS=['P42']
ANALYSIS_GROUPS = ['NR','DR','NRDR']
SUBGROUPS = ['annotation','classification']
SAMPLES = {
    "SC":{
        "P42":{
            "NR":{
                "NR1": os.path.join(DATA_DIR,"NR","AM3N207"),
                "NR2": os.path.join(DATA_DIR,"NR","AM3N208")
            },
            "DR":{
                "DR1": os.path.join(DATA_DIR,"DR","AM3N209"),
                "DR2": os.path.join(DATA_DIR,"DR","AM3N210")
            }
        }
    }
}

def init_paths(region='SC'):
    if not os.path.exists(REF_DATASETS):
        os.makedirs(REF_DATASETS)
    for timepoint in TIMEPOINTS:

        for group in ANALYSIS_GROUPS:

            for subgroup in SUBGROUPS:

                figure_dir = os.path.join(FIGURE_DIR,region,timepoint,group,subgroup)
                if not os.path.exists(figure_dir):
                    os.makedirs(figure_dir)
        
            data_dir = os.path.join(DATA_DIR,region,timepoint,group)
            if not os.path.exists(data_dir):
                os.makedirs(data_dir)

        model_dir = os.path.join(MODEL_DIR,region,timepoint)
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)
            print(model_dir)
        


                

def make_path(base='data',region='SC',timepoint="P42",group='NR',subgroup='annotation'):
    assert region in REGIONS and timepoint in TIMEPOINTS and group in ANALYSIS_GROUPS

    if base =="data":
        base_dir = DATA_DIR
        subgroup=""

    elif base=="figures":
        base_dir = FIGURE_DIR
        assert subgroup in SUBGROUPS

    else:
        raise AssertionError(f"Base {base} not a valid base directory")
    p = os.path.join(base_dir,region,timepoint, group,subgroup)
    assert os.path.exists(p)
    return p