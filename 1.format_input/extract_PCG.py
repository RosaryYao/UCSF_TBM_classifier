#!/usr/bin/python

import pandas as pd
import sys

inFile = sys.argv[1]
outFile = sys.argv[2]

#Input the .txt file you wish to run through the classifier#
df1 = pd.read_csv(inFile, sep = '\t')

#Input the ENSG_GeneSymbol_list.csv file that contains information on which genes are protein-coding#
df2 = pd.read_csv('data_files/ENSG_GeneSymbol_list.csv')


#Merge and extract the formatted file#

df3 = df1.merge(df2, how='inner', left_on='geneID', right_on='original')
df4 = df3[df3['gene_biotype'] == 'protein_coding']

formatted_data = df4.reset_index().drop('index',axis=1).drop('geneID',axis=1).set_index('ensembl').iloc[:,:-7]

formatted_data.to_csv(outFile)
