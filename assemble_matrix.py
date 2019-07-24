from __future__ import division
from sklearn.externals import joblib
from scipy import sparse
import numpy as np
import pandas as pd

precomputed_bins = joblib.load('./motifmap_cutoff_OpenChromChipSeq_protein_coding_data/huvec_precomputed_bins2')
mp_handler_result = joblib.load('./motifmap_cutoff_OpenChromChipSeq_protein_coding_NOHIC_data/mp_handler_result')
mp_handler_result = [(gene, sp_vector) for gene, sp_vector in mp_handler_result if not isinstance(sp_vector, str)]

q = pd.read_csv('./huvec_rnaseq_analysis/DE_genes_ALLGENES.csv')
z = q.loc[q.gene_biotype == 'protein_coding', 'ensembl_gene_id'] 
z = [i.lower() for i in z]

ml, tl = [], []

for i, j in enumerate(mp_handler_result):
    if i%100 == 0 or i == len(mp_handler_result)-1:
        if len(tl) != 0:
            ml.append(tl)
        tl = []
    tl.append(j)

batched_mp_handler = []
batched_genes = []

count = 0
for i, batch in enumerate(ml):    
   
    r = sparse.csr_matrix(np.ones(precomputed_bins.shape[1] - 1)) #seed final matrix with a row of ones; actual size, since binnums cut out
    gene_index_temp = []
    for gene, sp_vector in batch:
        count += 1
        print count
        if sp_vector == None or gene not in z:
            continue
        row = sp_vector.todense()
        r = sparse.vstack([r, row])
        gene_index_temp.append(gene)

    r = r.tocsc()[1:] #trim seed ones row
    batched_mp_handler.append(r)
    batched_genes.append(gene_index_temp)


r = sparse.csr_matrix(np.ones(precomputed_bins.shape[1] - 1)) #seed final matrix with a row of ones; actual size, since binnums cut out
final_gene_list = [] 
count = 0
print '***************'
for gene_batch, matrix_batch in zip(batched_genes, batched_mp_handler):
    count += 1
    print count
    r = sparse.vstack([r, matrix_batch])
    final_gene_list += gene_batch

print ('done')
r = r.tocsc()[1:] #trim seed ones row
print ('dumping')

joblib.dump(r, './motifmap_cutoff_OpenChromChipSeq_protein_coding_NOHIC_data/huvec_final_matrix_protein_coding')
joblib.dump(final_gene_list, './motifmap_cutoff_OpenChromChipSeq_protein_coding_NOHIC_data/huvec_final_matrix_genelist_protein_coding')
    

#once above line complete, vertically stack the vectors
#gene_index = []
#for gene, sp_vector in result:
#    if sp_vector == None:
#        continue
#    print gene, sp_vector.shape
#    #before stacking, densify
#    row = sp_vector.todense()
#    r = sparse.vstack([r, row])
#    gene_index.append(gene)
#
#r = r.tocsc()[1:] #trim seed ones row
#assert r.shape[0] == len(gene_index)
#return r, gene_index

