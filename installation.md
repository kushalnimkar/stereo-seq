Instructions taken from Berkeley CS 285 Fall 2023

1. Install conda, if you don't already have it, by following the instructions at [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

This install will modify the `PATH` variable in your bashrc.
You need to open a new terminal for that path change to take place (to be able to find 'conda' in the next step).

2. Create a conda environment that will contain python 3:
```
conda create -n stereo-seq python=3.10.4
```

3. activate the environment (do this every time you open a new terminal and want to run code):
```
source activate stereo-seq
```

4. Install the requirements into this conda environment
```
pip install -r requirements.txt
```

5. Allow your code to be able to see the 'stereo_seq' package installed by setup.py
```
cd <path_to_this directory>
pip install -e .
```
