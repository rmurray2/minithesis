from sklearn.externals import joblib
import collections
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#md = {'AUROC':[],
#      'Dataset':[],
#      'Set':[]}
#
#rand_labels = joblib.load('./HUVEC_up_AUROC_scores_shuffledy')
#true = joblib.load('./HUVEC_up_AUROC_scores')
#huvec_up = {'Random Labels':rand_labels, 'True Labels':true}
#for dataset, auroc_list in huvec_up.iteritems():
#    for auroc_score in auroc_list:
#        md['AUROC'].append(auroc_score)
#        md['Dataset'].append(dataset)
#        md['Set'].append('Up vs Rest')
#
#huvec_dn_rand = joblib.load('./HUVEC_dn_AUROC_scores_shuffledy')
#huvec_dn = joblib.load('./HUVEC_dn_AUROC_scores')
#huvec_dn_d = {'Random Labels':huvec_dn_rand, 'True Labels':huvec_dn}
#for dataset, auroc_list in huvec_dn_d.iteritems():
#    for auroc_score in auroc_list:
#        md['AUROC'].append(auroc_score)
#        md['Dataset'].append(dataset)
#        md['Set'].append('Down vs Rest')
#
#df = pd.DataFrame.from_dict(md)
#print df
#
#fig, axes = plt.subplots(1, 1)
#ax = sns.boxplot(x="Set", y="AUROC", hue="Dataset",data=df, ax=axes, width=.5, palette=['#597DBF', '#75BF71'], hue_order=['True Labels', "Random Labels"])
##ax = sns.catplot(x="Prior Threshold", y="AUROC", hue="Dataset", col="HiC", data=df, kind='box', height=4, aspect=2)
#ax.axhline(0.5, ls='--')
#plt.title('HUVEC siEPAS1 AUROCs (4-fold CV)')
#
#fig.savefig('HUVEC.png')

md = {'AUROC':[],
      'Dataset':[],
      'Set':[]}

rand_labels = joblib.load('./mef_up_AUROC_scores_shuffledy')
true = joblib.load('./mef_up_AUROC_scores')
mef_up = {'Random Labels':rand_labels, 'True Labels':true}
for dataset, auroc_list in mef_up.iteritems():
    for auroc_score in auroc_list:
        md['AUROC'].append(auroc_score)
        md['Dataset'].append(dataset)
        md['Set'].append('Up vs Rest')

mef_dn_rand = joblib.load('./mef_dn_AUROC_scores_shuffledy')
mef_dn = joblib.load('./mef_dn_AUROC_scores')
mef_dn_d = {'Random Labels':mef_dn_rand, 'True Labels':mef_dn}
for dataset, auroc_list in mef_dn_d.iteritems():
    for auroc_score in auroc_list:
        md['AUROC'].append(auroc_score)
        md['Dataset'].append(dataset)
        md['Set'].append('Down vs Rest')

df = pd.DataFrame.from_dict(md)
print df

fig, axes = plt.subplots(1, 1)
ax = sns.boxplot(x="Set", y="AUROC", hue="Dataset",data=df, ax=axes, width=.5, palette=['#597DBF', '#75BF71'], hue_order=['True Labels', "Random Labels"])
#ax = sns.catplot(x="Prior Threshold", y="AUROC", hue="Dataset", col="HiC", data=df, kind='box', height=4, aspect=2)
ax.axhline(0.5, ls='--')
plt.title('MEF shLKB1 (4-fold CV)')

fig.savefig('MEF.png')
#
