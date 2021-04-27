

def genes_to_xls(adata=None, filename='genes_hvg_regressed.xlsx', key='rank_genes_groups'):
    import xlsxwriter
    import pandas as pd
    filename = filename.replace('/','_')
    writer = pd.ExcelWriter(filename, engine='xlsxwriter')
    for i in adata.uns[key]['names'].dtype.names:
        d = {'gene':adata.uns[key]['names'][i],
             'score':adata.uns[key]['scores'][i],
             'mean':adata.var['mean'][adata.uns[key]['names'][i]]} #note the mean comes with an index.
        df = pd.DataFrame(data=d)
        i=i.replace('/', '_')
        df.to_excel(writer,sheet_name=i)
    writer.save()
    return
