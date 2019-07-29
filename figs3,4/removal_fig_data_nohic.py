from __future__ import division
import pandas as pd
import operator
from sklearn import cross_validation
import numpy as np
from sklearn.externals import joblib
from sklearn import feature_selection
import xgboost
from sklearn import preprocessing
from sklearn.metrics import roc_auc_score, average_precision_score, precision_score, recall_score

name_to_enst = {}
enst_to_ensg = {}

with open('./ensemblToGeneName.txt') as infile:
    for line in infile:
        sline = line.split('\t')
        gene_name = sline[1].lower().strip()
        trxpt_name = sline[0]
        if gene_name not in name_to_enst:
            name_to_enst[gene_name] = [trxpt_name]
        else:
            name_to_enst[gene_name].append(trxpt_name)

with open('./ensGene.txt') as infile:                                                                                                                                                                                                                           
    for line in infile:                                                                                                                                                                                                                                                
        sline = line.split('\t')                                                                                                                                                                                                                                                
        if sline == []: break                                                                                                                                                                                                                                                   
        gene_name = sline[12].lower().strip()                                                                                                                                                                                                                                   
        trxpt_name = sline[1]                                                                                                                                                                                                                                           
        if trxpt_name not in enst_to_ensg:                                                                                                                                                                                                                                        
            enst_to_ensg[trxpt_name] = [gene_name]                                                                                                                                                                                                                               
        else:                                                                                                                                                                                                                                                                   
            enst_to_ensg[trxpt_name].append(gene_name)                                                                                                                                                                                                                           
    
name_to_ensg = {k:enst_to_ensg[v[0]][0] for k,v in name_to_enst.iteritems()}
ensg_to_name = {v:k for k,v in name_to_ensg.iteritems()}

for i, j in name_to_ensg.iteritems():
    print i, j
    break
  
for i, j in ensg_to_name.iteritems():
    print i, j
    break  


pnas_df = pd.read_csv('./pnas.1802973116.sd02.txt', delimiter='\t')
pnas_d = {j.Gene_Name.lower():j.DE_Prior_Rank for i, j in pnas_df.iterrows()}

cutoffs = [1, .9, .8, .7, .6, .5, .4, .3, .2, .1]

y = joblib.load('./motifmap_cutoff_OpenChromChipSeq_protein_coding_NOHIC_data/y/up_AllMags_.05padj')
genes = joblib.load('./motifmap_cutoff_OpenChromChipSeq_protein_coding_NOHIC_data/huvec_final_matrix_genelist_protein_coding')
X = joblib.load('./motifmap_cutoff_OpenChromChipSeq_protein_coding_NOHIC_data/huvec_final_matrix_protein_coding')
# feature_labels = joblib.load('./motifmap_cutoff_OpenChromChipSeq_protein_coding_data/huvec_feature_labels')
selector = feature_selection.VarianceThreshold(.8 * (1 - .8))
X_new = selector.fit_transform(X)
print X_new.shape

#     feature_labels_new = np.array(feature_labels)
#     feature_labels_new = feature_labels_new[selector.get_support() == True]
#     X_new = X

scale = preprocessing.StandardScaler(with_mean=False)
X_new = scale.fit_transform(X_new)
X_new = X_new.todense()
X_new = np.array(X_new)
    
mr = {i:[] for i in cutoffs}

for cutoff in cutoffs:
    print ("")
    print 'cutoff', cutoff
    new_y = []
    for gene, label in zip(genes, y):
#         print gene
        if gene not in ensg_to_name:
            new_y.append(label)
            continue
        if ensg_to_name[gene] not in pnas_d:
            new_y.append(label)
            continue
        
        #if the prior of the gene is over the cutoff, set label to zero
        if pnas_d[ensg_to_name[gene]] >= cutoff:
            new_y.append(0)
            
        else:
            new_y.append(label)
        
    new_y = np.array(new_y)
    print np.bincount(new_y)
    
    count = 0
    while count < 25:
        clf = xgboost.XGBClassifier(n_estimators=100, n_jobs=12)

        skf = cross_validation.StratifiedKFold(new_y, n_folds=4, shuffle=True)
        for i, (train, test) in enumerate(skf):
            preds = clf.fit(X_new[train], new_y[train], 
                            eval_metric='auc').predict_proba(X_new[test])
            auroc = roc_auc_score(new_y[test], preds[:,1])
            print auroc
            mr[cutoff].append(auroc)
        count += 1
        
joblib.dump(mr, './HUVEC_AUROC_removal_titration_results_nohic')         
