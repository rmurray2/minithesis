from sklearn.externals import joblib
import collections
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#rand_probs = joblib.load('./HUVEC_AUROC_removal_titration_results_randomprobs')
#real = joblib.load('./HUVEC_AUROC_removal_titration_results')
#rand_labels = joblib.load('./HUVEC_AUROC_removal_titration_results_randomlabels')

rand_probs_nohic = joblib.load('./HUVEC_AUROC_removal_titration_results_randomprobs_nohic')
real_nohic = joblib.load('./HUVEC_AUROC_removal_titration_results_nohic')
rand_labels_nohic = joblib.load('./HUVEC_AUROC_removal_titration_results_randomlabels_nohic')

md = {'AUROC':[],
      'Prior Threshold':[],
      'Dataset':[],
      'HiC':[]}

#ind = {'Random Priors':rand_probs, 'Random Labels':rand_labels,'True Priors + Labels':real}
ind_nohic = {'Random Priors':rand_probs_nohic, 'Random Labels':rand_labels_nohic,'True Priors + Labels':real_nohic}
#
#for dataset, data in ind.iteritems():
#        for cutoff, auroc_list in data.iteritems():
#            for auroc_score in auroc_list:
#                md['AUROC'].append(auroc_score)
#                md['Prior Threshold'].append(cutoff)
#                md['Dataset'].append(dataset)
#                md['HiC'].append('+')
#
for dataset, data in ind_nohic.iteritems():
        for cutoff, auroc_list in data.iteritems():
            for auroc_score in auroc_list:
                md['AUROC'].append(auroc_score)
                md['Prior Threshold'].append(cutoff)
                md['Dataset'].append(dataset)
                md['HiC'].append('-')

df = pd.DataFrame.from_dict(md)

fig, axes = plt.subplots(1, 1)
ax = sns.boxplot(x="Prior Threshold", y="AUROC", hue="Dataset",data=df, ax=axes, width=.5, order=[1, .9, .8, .7, .6, .5, .4, .3, .2, .1], palette='muted', hue_order=['True Priors + Labels', 'Random Priors', "Random Labels"])
#ax = sns.catplot(x="Prior Threshold", y="AUROC", hue="Dataset", col="HiC", data=df, kind='box', height=4, aspect=2)
ax.axhline(0.5, ls='--')
plt.title('AUROC as a function of Prior Threshold: HUVEC EPAS1 KD (-HiC)\nUp Vs Rest')

fig.savefig('AURUC_Prior-HiC.png')

