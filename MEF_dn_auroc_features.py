import pandas as pd
import operator
from sklearn import cross_validation
import numpy as np
from sklearn.externals import joblib
from sklearn import feature_selection
import xgboost
from sklearn import preprocessing
from sklearn.metrics import roc_auc_score
import numpy as np

X = joblib.load('./mef_data/mef_final_matrix')
feature_labels = joblib.load('./mef_data/mef_feature_labels')
gene_labels = joblib.load('./mef_data/mef_final_matrix_genes') 
dfo = pd.read_table('./mef_data/gene_exp.diff')
print dfo.columns
y = []

c = 0
up_genes = []
for i, row in dfo.iterrows():
    if (row.status == 'OK') and  \
    (row['log2(fold_change)'] <= .75) and  \
    (row.p_value <= .05) and \
    row.gene.lower() in gene_labels:
          up_genes.append(row.gene.lower())
          c += 1
print c

for gene in gene_labels:
    if gene in up_genes:
        y.append(1)
    else:
        y.append(0)

################################

#feature selection
selector = feature_selection.VarianceThreshold(.8 * (1 - .8))
X_new = selector.fit_transform(X)
print X_new.shape

#scaling
scale = preprocessing.StandardScaler(with_mean=False)
X_new = scale.fit_transform(X_new)
X_new = X_new.todense()
X_new = np.array(X_new)

y = np.array(y)
print np.bincount(y)

#####CV real data

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
    
joblib.dump(aurocs_real, './mef_dn_AUROC_scores')

#train on all features, scaled. use probability-shifted y values

scale = preprocessing.StandardScaler(with_mean=False)
X_new = scale.fit_transform(X)
X_new = X_new.todense()
X_new = np.array(X_new)

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

joblib.dump(od, './mef_dn_features')

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
    
joblib.dump(aurocs_random, './mef_dn_AUROC_scores_shuffledy')
