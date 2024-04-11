Shared private github repo for Single-cell and stereo-seq data analysis of mouse in control and dark reared conditions

# FAQ

## Error messages:
"Missing classes: {missing_classes}, missing_clusters = {missing_clusters}", when running 
```create_subtypes(adata,groupby,subgroups,key='leiden',key_added='class_subtypes',output_f="")```
    This error means that there are classes that have been found by leiden but isn't in the list under the parameter subgroups.