# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 23:17:53 2019

@author: Лена
"""

#from intermine.webservice import Service
#import pandas as pd
#def thalemine_query(gene):
#    service = Service("https://apps.araport.org:443/thalemine/service")
#    query = service.new_query("Gene")
#    query.add_view(
#            "primaryIdentifier", "RNASeqExpressions.expressionLevel",
#            "RNASeqExpressions.experiment.SRAaccession",
#            "RNASeqExpressions.experiment.tissue", "RNASeqExpressions.unit")
#    query.add_sort_order("Gene.RNASeqExpressions.expressionLevel", "DESC")
#    query.add_constraint("Gene", "LOOKUP", gene, "A. thaliana", code = "A")
#    column_name = ["gene_id", "expressionLevel", "SRAaccession", "experiment_tissue", "unit"]
#    temp_df = []
#    for row in query.rows():
#        raw_row = []
#        raw_row.append(row["primaryIdentifier"])
#        raw_row.append(row["RNASeqExpressions.expressionLevel"])
#        raw_row.append(row["RNASeqExpressions.experiment.SRAaccession"])
#        raw_row.append(row["RNASeqExpressions.experiment.tissue"])
#        raw_row.append(row["RNASeqExpressions.unit"]) 
#        temp_df.append(raw_row)
##    return(temp_df)
#    data = pd.DataFrame(temp_df, columns = column_name)
#    data.to_csv(gene + '_RNAseq' + '.csv', index = False)
#    print(gene, 'data loaded succesfully')
#    return data
##a = thalemine_query('AT1G05010')
#gene = pd.read_fwf('gene_tf.txt', header=None, names = ['gene_id', 'trivial'])
#i = 0
#for gene_id in gene.loc[:, 'gene_id']:
#    thalemine_query(gene_id)
#print('all right!')
##gene.columns = ['gene_id', 'trivial']
##row["RNASeqExpressions.expressionLevel"], \
##        row["RNASeqExpressions.experiment.SRAaccession"], \
##        row["RNASeqExpressions.experiment.tissue"], row["RNASeqExpressions.unit"]
a = range(1392)
for i in a[::232]:
    print(i)