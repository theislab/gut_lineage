import pandas as pd
import patsy
import sys
import numpy.linalg as la
import numpy as np
import scipy as sci
from matplotlib import rcParams
import matplotlib.pyplot as pl

def round_sig(x, sig=2, small_value=1.0e-100):
    expo = int(np.floor(np.log10(max(abs(x), abs(small_value)))))
    base = round(x * 10**-expo, sig)
        
    return(base, expo)

def bar_frequency(adata=None,cluster='louvain', group=None,celltype1=None,celltype2=None, do_test=False):
    


    sub_cells = np.in1d(adata.obs[group], [celltype1])
    print(celltype1)
    print(np.sum(sub_cells))
    adata_filt=adata[sub_cells,:]
    data = adata_filt.obs[cluster]
    group_names=np.unique(adata.obs[cluster])
    counts_control = pd.value_counts(data,normalize=True)
    counts_control = counts_control[group_names]
    counts_control[np.isnan(counts_control)] = 0
    counts_control2 = pd.value_counts(data,normalize=False)
    counts_control2 = counts_control2[group_names]
    counts_control2[np.isnan(counts_control2)] = 0
    
    n_groups = len(group_names)

    sub_cells = np.in1d(adata.obs[group], [celltype2])
    print(celltype2)
    print(np.sum(sub_cells))
    adata_filt=adata[sub_cells,:]
    data = adata_filt.obs[cluster]
    counts_enriched = pd.value_counts(data,normalize=True)
    counts_enriched = counts_enriched[group_names]
    counts_enriched[np.isnan(counts_enriched)] = 0
    counts_enriched2 = pd.value_counts(data,normalize=False)
    counts_enriched2 = counts_enriched2[group_names]
    counts_enriched2[np.isnan(counts_enriched2)] = 0
    
    #compute p-value for each cluster with Fisher's exact test
    p_val = np.zeros(n_groups)
    for x in range(0, n_groups):
        non_c_enriched = counts_enriched2[np.arange(n_groups)!=x].sum()
        non_c_control = counts_control2[np.arange(n_groups)!=x].sum()
        testmat = np.array([[counts_control2[x], non_c_control],[counts_enriched2[x], non_c_enriched]], dtype=int)
        p_val[x]=min(1, sci.stats.fisher_exact(testmat, alternative='two-sided')[1]* n_groups) #Bonferroni correction
    if do_test:
        res_order = np.argsort(p_val)

    else:
        res_order =np.arange(n_groups)
    #print(counts_control[res_order])
    #print(counts_enriched[res_order])
    #print(res_order)
    #print(p_val)
    
    rcParams['text.usetex']=True
    fig, ax = pl.subplots()

    index = np.arange(n_groups)
    bar_width = 0.35

    opacity = 0.6
    error_config = {'ecolor': '0.3'}
    
    if type(celltype1)==str:
        label1 = celltype1.replace('_', ' ')

    else:
        label1='Control'
    
    if type(celltype2)==str:
        label2 = celltype2.replace('_', ' ')
    else:
        label2='Enriched'
    
    groupp = group.replace('_', ' ')
    rects1 = pl.bar(index, counts_control[res_order], bar_width,
                     alpha=opacity,
                     color='grey',
                     error_kw=error_config,
                     label=label1)

    rects2 = pl.bar(index + bar_width, counts_enriched[res_order], bar_width,
                     alpha=opacity,
                     color='blue',
                     error_kw=error_config,
                     label=label2)
    
    
    for i in index:
        ord_i = res_order[i]
        rounded_p_val = round_sig(p_val[ord_i], sig=2)
        if rounded_p_val[1]==0:
            str_p_val ='$1.0$'
        else:        
            str_p_val = '$' + str(rounded_p_val[0]) + '\cdot 10^{' + str(rounded_p_val[1]) + '}$'
        
        pl.text(i+0.5*bar_width, max(counts_control[ord_i], counts_enriched[ord_i]), str_p_val, 
                horizontalalignment='center',
                verticalalignment='center')
        
    pl.xlabel('Cluster')
    pl.title(groupp)
    pl.ylabel('\% Cells')
    pl.xticks(index + bar_width / 2, group_names[res_order], rotation='vertical')
    
    #pl.title('Scores by group and gender')
    pl.tight_layout()
    pl.legend()
    pl.show()
    rcParams['text.usetex']=False #disable LaTeX usage after plotting as it messes with other plot functions
    return()
