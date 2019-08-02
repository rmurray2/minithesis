from __future__ import division
import operator
import pandas as pd
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

y = joblib.load('./motifmap_cutoff_OpenChromChipSeq_protein_coding_data/y/dn_AllMags_.05padj')
genes = joblib.load('./motifmap_cutoff_OpenChromChipSeq_protein_coding_data/huvec_final_matrix_genelist_protein_coding')
X = joblib.load('./motifmap_cutoff_OpenChromChipSeq_protein_coding_data/huvec_final_matrix_protein_coding')
feature_labels = joblib.load('./motifmap_cutoff_OpenChromChipSeq_protein_coding_data/huvec_feature_labels')

#scaling
scale = preprocessing.StandardScaler(with_mean=False)
X_new = scale.fit_transform(X)
X_new = X_new.todense()
X_new = np.array(X_new)

y = np.array(y)
print np.bincount(y)

aurocs_real = []
count = 0
while count < 50:

    skf = cross_validation.StratifiedKFold(y, n_folds=4, shuffle=True)
    for i, (train, test) in enumerate(skf):
        clf = xgboost.XGBClassifier(n_estimators=100, n_jobs=12)
        preds = clf.fit(X_new[train], y[train], 
                        eval_metric='auc').predict_proba(X_new[test])
        auroc = roc_auc_score(y[test], preds[:,1])
        print auroc
        aurocs_real.append(auroc)
    count += 1
    
joblib.dump(aurocs_real, './HUVEC_dn_AUROC_scores')

#train on all features, scaled. use probability-shifted y values

clf = xgboost.XGBClassifier(n_estimators=200,  n_jobs=12)
clf.fit(X_new, y, eval_metric='auc')       

res = clf.get_booster().get_score(importance_type='total_gain')
importance = sorted(res.items(), key=operator.itemgetter(1))

od = {}
for i, j in importance:
    idx = int(i[1:])
    name = feature_labels[idx]
    print name, j
    od[name] = j

joblib.dump(od, './HUVEC_dn_features')

########CV shuffled DATA


aurocs_random = []
count = 0
while count < 50:

    np.random.shuffle(y)
    skf = cross_validation.StratifiedKFold(y, n_folds=4, shuffle=True)
    for i, (train, test) in enumerate(skf):
        clf = xgboost.XGBClassifier(n_estimators=100, n_jobs=12)
        preds = clf.fit(X_new[train], y[train], 
                        eval_metric='auc').predict_proba(X_new[test])
        auroc = roc_auc_score(y[test], preds[:,1])
        print auroc
        aurocs_random.append(auroc)
    count += 1
    
joblib.dump(aurocs_random, './HUVEC_dn_AUROC_scores_shuffledy')
