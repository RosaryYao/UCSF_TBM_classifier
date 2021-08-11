#!/usr/bin/python

import pandas as pd
import numpy as np
import random
import sys
import csv
import os
from sklearn.svm import LinearSVC
from sklearn import svm
from sklearn.model_selection import GridSearchCV
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
import sklearn.metrics as metrics
import warnings

warnings.filterwarnings("ignore")


####User input the files#### 

inFile = sys.argv[1]
outFile = sys.argv[2]

dataframe1 = pd.read_csv('data_files/training_validation.csv').set_index('ensembl') 
dfb = pd.read_csv('data_files/metadata.csv').set_index('samples')

df_TEST = pd.read_csv(inFile).set_index('ensembl')

####Functions for MLC######

def pick_ref_sample(x):
    
    """Pick a reference sample to calculate normalization factors for input sample. Output from this function is fed as input into calcNormFactors_py()"""
    
    #Remove genes that are zeros throughout
    num = x.shape[1]
    headers = x.columns
    x['sum'] = x.sum(axis=1)
    y = x[~x['sum'].isin([0])]
    z = y.iloc[:,0:num].T
    xx = z 
    
    #Get cpm
    xx['sum'] = xx.sum(axis=1)
    zz = (xx.loc[:, xx.columns != "sum"].div(xx["sum"], axis=0).T).iloc[:,:(num)]
    
    #Get 75% quant table
    quant_75_table = zz.quantile([.75], axis = 0).T
    average_quant_75_table = float(zz.quantile([.75], axis = 0).T.sum()/len(zz.quantile([.75], axis = 0).T))
    
    #Pick reference sample using edgeR method
    
    Ref_sample_table = abs(quant_75_table-average_quant_75_table).sort_values(by=0.75)
    Ref_sample = Ref_sample_table.index[0]
    
    
    return (zz, Ref_sample)


def calcNormFactors_py(zz, Ref_sample, output_scaling_factor_file,dataframe_OG):
    
    """Python implementation of edgeR's calcNormFactors."""
    
    #Calculate log(2) ratio's of reference sample/input sample
    Ratio_zz = 1/(zz.loc[:].div(zz[Ref_sample], axis=0))
    mat1 = np.log2(Ratio_zz)
    
    #Array of sample IDs
    all_samples = mat1.columns
    
    textfile = open(output_scaling_factor_file,"a")

    for i in range(len(all_samples)): 
        if all_samples[i] == Ref_sample:
            #scaling factor for Ref_sample is 1
            print (all_samples[i], 1, file=textfile,sep=",")
        else:
            #generate 2 matrices 1) Genes in ref and input samples have finite values, 2) Geometric mean of genes in ref and input samples
            mat_t = mat1[[all_samples[i], Ref_sample]]
            mat_t1 = mat_t[(~mat_t[all_samples[i]].isin(['NaN', 'inf','-inf'])) & (~mat_t[Ref_sample].isin(['NaN', 'inf','-inf']))]
            mat_t1_sort = mat_t1.sort_values(by=all_samples[i])
    
            t1 = ((mat1[all_samples[i]] + mat1[ Ref_sample] )/ 2).to_frame()
            mat_t2 = t1[(~t1[0].isin(['NaN', 'inf','-inf',0]))] #in addition to filtering out genes that are infinity, remove genes that are 0 as well
            mat_t2_sort = mat_t2.sort_values(by=0)

            #filter out top and bottow 30% genes of mat_t1_sort, and top and bottom 5% genes of mat_t2_sort
            num30 = int(np.round(len(mat_t1_sort)*0.30))
            num05 = int(np.round(len(mat_t2_sort)*0.05))

            keep_after30 = mat_t1_sort.index[num30:(len(mat_t1_sort.index))-num30]
            keep_after05 = mat_t2_sort.index[num05:(len(mat_t2_sort.index))-num05]
            
            #Keep genes common to both keep_after30 and keep_after05 lists
            genes_common = np.intersect1d(keep_after30,keep_after05)
                        
                
            #If the intersection between them is zero, drop until it is not
            if len(genes_common) == 0:
                j = 30
                while len(genes_common) == 0:
                    num30 = int(np.round(len(mat_t1_sort)*j)/100)
                    num05 = int(np.round(len(mat_t2_sort)*0.05))

                    keep_after30 = mat_t1_sort.index[num30:(len(mat_t1_sort.index))-num30]
                    keep_after05 = mat_t2_sort.index[num05:(len(mat_t2_sort.index))-num05]
                    genes_common = np.intersect1d(keep_after30,keep_after05) 
                    j = j - 1
                    if j == 1:
                        break
                
            #calculate scaling factor
            log_ratio_sample = mat1[all_samples[i]].loc[genes_common].to_frame()
            original_sample = dataframe_OG[all_samples[i]].loc[genes_common].to_frame()
           
            if len(log_ratio_sample) == 0:
                scaling_factor = 0
            else:
                weighted_avg = float((log_ratio_sample * original_sample).sum())/float(original_sample.sum())
                scaling_factor = 2**(-weighted_avg)
            print(all_samples[i], scaling_factor, file=textfile,sep=",")
        
    textfile.close()
    return 


def normalize_matrix_scaling_logcpm(df,scaling_factors,prior_count):
    
    """Peform the log CPM transformation using scaling factors from edgeRs method"""
    #Read in input raw count dataframe, and scaling factor CSV file obtained using edgeR method
    test = df.T.reset_index()
    scaling = pd.read_csv(scaling_factors, header = None)
    scaling.drop_duplicates(keep='first', inplace=True)
    test_x = test.merge(scaling,how='left', left_on='index', right_on=0).drop(0,axis=1).set_index('index')
    
    #scale matrix with normlaization factors
    new_scale = test_x.loc[:, test_x.columns != 1].div(test_x[1], axis=0)
    
    #CPM normalization
    new_scale['sum'] = test_x.sum(axis=1)
    scaled = (new_scale.loc[:, new_scale.columns != "sum"].div(new_scale["sum"], axis=0))*1000000
    
    #log2 transform using user inputted prior count
    value_for_CF = np.log2(scaled+prior_count).T

    return value_for_CF


####MLC predictions######

f = open(outFile, "w")
writer = csv.DictWriter(
    f, fieldnames=["Sample ID", "Predicted value"])
writer.writeheader()
f.close()

textfile = open(outFile,"a")

BS = 100 #How many times to bootstrap

for i in np.arange(BS):
    
    ####split the dataset into 80% to test your classifier on the current classifier###

    n = dataframe1.shape[1]
    perm = np.random.permutation(n)
    a = int(np.round(0.8 * n))

    new_perm = dataframe1.iloc[:,perm]
    train_x = new_perm.iloc[:,:a]

    
    ####perform normalization####
    
    #pick the reference samples from the training data#
    norm_trainx, ref = pick_ref_sample(train_x)
    
    #calculate scaling factors for the training data#
    calcNormFactors_py(norm_trainx, ref, 'scaling_factors.csv', dataframe1)
    
    #convert validation matrix to a format suitable for calcNormFactors_py function#
    test_TrainGenes_x = norm_trainx[ref].to_frame().merge(df_TEST, how='left', left_index = True, right_index = True).drop(ref, axis = 1)
    test_TrainGenes = test_TrainGenes_x.merge(new_perm[ref].to_frame(), how='left', left_index = True, right_index = True).T
    
    test_TrainGenes['sum'] = test_TrainGenes.sum(axis=1)
    norm_testx = (test_TrainGenes.loc[:, test_TrainGenes.columns != "sum"].div(test_TrainGenes["sum"], axis=0).T)
    
    #calculate scaling factors for the training data#
    calcNormFactors_py(norm_testx, ref, 'scaling_factors.csv', df_TEST)
    
    #calculate normalized training and validation matices#
    new_perm_x = normalize_matrix_scaling_logcpm(new_perm,'scaling_factors.csv',75)
    train_norm = new_perm_x.iloc[:,:a].dropna(axis=1)
    test_norm = normalize_matrix_scaling_logcpm(df_TEST,'scaling_factors.csv',75).dropna(axis=1)
    
       
    #remove scaling_factors.csv file
    os.remove("scaling_factors.csv")
    
    
    ####feature selection####
    
    #merge the key with the normalized training dataset
    train_norm_merged = train_norm.T.merge(dfb, how='left', left_index = True, right_index = True)
    train_norm_fit = train_norm_merged.iloc[:,0:-1]
    train_norm_fit_key = train_norm_merged.iloc[:,-1]

    #perform the chi-square analysis
    bestfeatures = SelectKBest(score_func=chi2, k=40)
    fit = bestfeatures.fit(train_norm_fit,train_norm_fit_key)
    dfscores = pd.DataFrame(fit.scores_)
    dfcolumns = pd.DataFrame(train_norm_fit.columns)
    
    #concat genes and scores into one dataframe, and pick top 500 features
    featureScores = pd.concat([dfcolumns,dfscores],axis=1)
    featureScores.columns = ['Genes','Score']
    s_featureScores = featureScores.sort_values(['Score'], ascending=False)
    univariate_features = s_featureScores[:40]['Genes'].values
    
    #generate the new training and validation dataset based on univariate_features
    train_norm_x = train_norm[train_norm.index.isin(univariate_features)]
    train_norm_x_merged = train_norm_x.T.merge(dfb, how='left', left_index = True, right_index = True)
    
    #generate the final training and test sets
    TB_doub = train_norm_x_merged[train_norm_x_merged['key']==1]
    training1 = train_norm_x_merged.drop('key',axis=1)
    training2 = TB_doub.drop('key',axis=1) 

    training = training1.T.merge(training2.T, how='left', left_index=True, right_index=True).T

    Y1 = train_norm_x_merged.loc[:,'key']
    Y2 = TB_doub.loc[:,'key']

    Y = np.concatenate((Y1, Y2))
    
    test = test_norm[test_norm.index.isin(univariate_features)].T

    
    
    ####SVM####
    
    lsvc = LinearSVC(C=1/10, penalty="l1", dual=False).fit(training.values,Y)
    model = SelectFromModel(lsvc, prefit=True)
    predicted = lsvc.predict(test.values) #predicting on the remaining samples

    idx = test.index
    SVMpred = predicted
    
    for k in range(len(idx)):
        print (idx[k], SVMpred[k], file=textfile,sep=",")
        
textfile.close()
