from __future__ import division
from sklearn.preprocessing import MinMaxScaler
from sklearn.externals import joblib
import multiprocessing
from scipy import sparse
import pandas as pd
import numpy as np
import datetime
import bisect
import gc
import random

'''
given a cufflinks gene_exp.diff file, produce a sparse matrix
of hic-weighted tfbs/chromatin features from precomputed bins
'''
# #load
#     precomputed bin vectors
#get min/max frequencies
# #create
#     bin_mapping_dict
#     name_to_ens
#     ens_to_pos
#     bin number : precomputed index dict
#precomputed_bins = joblib.load('./huvec_precomputed/huvec_precomputed_bins')
precomputed_bins = joblib.load('./huvec_precomputed/huvec_precomputed_bins2')
bin_ids = precomputed_bins[:,0].todense()
bin_ids = np.array(bin_ids)
#not have any  predicted tfbs

#given bin number, look up row index in precomputed matrix
idx_dict = {}
for index, binnum in enumerate(bin_ids):
    idx_dict[binnum[0]] = index

freqs = []
#with open('./HiC/imr90_tnfa/imr90_tnfa_results/hic_results/matrix/imr90_tnfa/iced/30000/imr90_tnfa_30000_iced.matrix') as infile:
#with open('./HiC/mef_battulin/hic_results_30K/matrix/mef_exp/iced/30000/mef_exp_30000_iced.matrix') as infile:
with open('./huvec_out/hic_results/matrix/huvec/iced/30000/huvec_30000_iced.matrix', "rb") as infile:
    for line in infile:
        sline = line.split('\t')
        freq = float(sline[-1])
        freqs.append(freq)

max_freq, min_freq = max(freqs), min(freqs)
################
bin_mapping_dict = {}
name_to_ens = {} #{genename : [ENST1, ENST2]}
ens_to_pos = {} # {ENST : (chr1, 92349) }

#bin mapping dict (d) is of the form: { chr1 : [ starting bin for key chrom, last position value in key chrom ]

#binmapping_f = open('./HiC/huvec_validation/huvec_results/hic_results/matrix/huvec/raw/30000/huvec_30000_abs.bed').read().split('\n')
#binmapping_f = open('./HiC/imr90_tnfa/imr90_tnfa_results/hic_results/matrix/imr90_tnfa/raw/30000/imr90_tnfa_30000_abs.bed').read().split('\n')
#binmapping_f = open('./HiC/mef_battulin/hic_results_30K/matrix/mef_exp/raw/30000/mef_exp_30000_abs.bed').read().split('\n')
binmapping_f = open('./huvec_out/hic_results/matrix/huvec/raw/30000/huvec_30000_abs.bed').read().split('\n')

#create dict
for line in binmapping_f[:-1]:
    sline = line.split('\t')
    if sline[0] not in bin_mapping_dict:
        bin_mapping_dict[sline[0]] = [0, int(sline[3])]

    bin_mapping_dict[sline[0]] = [int(sline[2]), bin_mapping_dict[sline[0]][1]]

with open('./ensGene.txt') as infile:
    for line in infile:

        sline = line.split('\t')
        if sline == []: break

        gene_name = sline[12].lower().strip()
        trxpt_name = sline[1]
        if gene_name not in name_to_ens:
            name_to_ens[gene_name] = [trxpt_name]
        else:
            name_to_ens[gene_name].append(trxpt_name)

        strand = sline[3]
        ens = sline[1]
        chrom = sline[2]
        first_pos = int(sline[4])
        second_pos = int(sline[5])

        if strand == '+':
            ens_to_pos[ens] = (chrom, first_pos)

        if strand == '-':
            ens_to_pos[ens] = (chrom, second_pos)

print 'done initializing'
#pseudocode
# final_matrix = []
# for each row:
#     get the pvalue
#
#     #get the gene's bin number
#     selfbin = find_bin(gene) #gene is a gene name

#     #get all bins interacting with gene-bin
#     intxng_bins = [] #get neighbors from networkx object
#     intxng_bin_freq = [] #get neighbor interaction frqeuency from networkx object

#     #if no interacting bins, just use the self bin
#     if no interacting bins:
#         get self bin vector and append to final matrix
#         continue

#     #for each interacting bin, get the bins data vector
#     interacting_bin_vectors = []
#     for binnum in intxng_bins:
#         bv = get_bin_vector(binnum)
#         interacting_bin_vectors.append(bv.todense())

#     #combine vectors using frequency weights; remember ot sum freqs and divide by N
#     final_vector = combine_vectors(interacting_bin_vectors, intxng_bin_freq)
#     final_matrix.appenD(final_vector)

def get_gene_list(gene_diff):
    q = pd.read_csv(gene_diff)
    z = q.loc[q.gene_biotype == 'protein_coding', 'ensembl_gene_id'].tolist()
    genes = [i.lower() for i in z]

#    df = pd.read_csv(gene_diff)
#
#    genes = df['ensembl_gene_id'].tolist()
#    genes = [i.lower() for i in genes]

    final_list = []
    for i in genes:
        if i in name_to_ens:
            final_list.append(i) # for each transcript take the first name that maps to an enstr

    final_list = list(set(final_list)) #remove dupes
    return final_list

def findbin(chrom, start, binsize):
    '''
    given chr position, binsize and prepared chr dict, returns bin of position
    '''
    return bisect.bisect(range(0, bin_mapping_dict[chrom][0], binsize), start) + bin_mapping_dict[chrom][-1] - 1

def gene_to_bins(gene_name):
    '''
    takes gene name, lower&stripped and returns the bin number(s) as a list of ints
    '''

    gene_name = gene_name.lower()

    ensnames = name_to_ens[gene_name]
    positions = []
    for transcript in ensnames:
        pos = ens_to_pos[transcript]
        if pos[0] in bin_mapping_dict:
            positions.append(pos)

    binnums = []
    for pos in positions:
        binnum = findbin(pos[0], pos[1], 30000)
        binnums.append(binnum)

    #to get list of unique bins, uncomment next line
#    binnums = list(set(binnums))

    return binnums

def find_interacting_bins_str(binnum):
    '''
    STR comparison
    def to take in bin number, give back
    all interacting bins and freqs as
    {bin : freq}

    '''
    binnum_str = str(binnum)
    rd = {}
#    with open('./HiC/huvec_validation/huvec_results/hic_results/matrix/huvec/iced/30000/huvec_30000_iced.matrix', "rb") as infile:
#    with open('./HiC/imr90_tnfa/imr90_tnfa_results/hic_results/matrix/imr90_tnfa/iced/30000/imr90_tnfa_30000_iced.matrix') as infile:
    with open('./huvec_out/hic_results/matrix/huvec/iced/30000/huvec_30000_iced.matrix') as infile:
        for line in infile:
            if binnum_str in line: #only process line if the query binnum is somewhere on it

                sline = line.split('\t')
                j, k = int(sline[0]), int(sline[1])
                if j > binnum: #since first col sorted and second never < first, rest of file not worth iterating
                    return rd
                if binnum == j:
                    l = float(sline[2])
                    rd[k] = l

                elif binnum == k:
                    l = float(sline[2])
                    rd[j] = l
    return rd

def mp_worker(gene):
    '''
    main func that creates final data vector, row by row for each gene
    in is gene name
    out is sparse tfbs data vector for that gene
    the vector is a combination of precomputed data vectors corresp to hic
    interacting bins, weighted by intxn freq
    '''
    start =datetime.datetime.now()
    try:
        #convert gene name to list of bin numbers (coresp to diff trxpt TSSs)
        binnums = gene_to_bins(gene.lower())
        assert binnums != []

        #pick bin that includes the most transcripts
        #unfortunately cufflinks doesn't specify tss for rows, so have to
        #make a best guess
        chosen_bin = max(set(binnums), key=binnums.count)
        print binnums
        print chosen_bin

        #for the chosen bin, find all interacting bins
        #this step might be improved somehow
        #iced matrix file for each chosen bin
        inxtng_bins = find_interacting_bins_str(chosen_bin)

        #set freq of self bin just as high as the freq max
        #rationale  is that tfbs/chommods around promoter are
        #v imp features in detemrining diff exp
        if chosen_bin not in inxtng_bins: #if self not in inxtng bin list
            #add it, and set freq to max of rest of freqs
            #or 1 if no other interacting bins
            if len(inxtng_bins) == 0:
                mv = 1
            else:
                mv = max(inxtng_bins.values()) #TODO maybe multiple this by somethign?

            inxtng_bins[chosen_bin] = mv

        #get the bin data vectors for the interacting bins  and densify
        interacting_bin_vectors = []
        interacting_bin_freqs = []

        #skip interacting bins that have no precomputed data vector
        #the bin might not have any pred tfbs sites

        for binn, freq in inxtng_bins.iteritems():
            if binn not in idx_dict: #skip interacting bins that have no precomputed data vector
                continue
            index = idx_dict[binn] #get precomputed bin matrix index (row)
            bin_data = precomputed_bins.getrow(index)
            bin_data = np.array(bin_data.todense()) #densify
            bin_data = bin_data[0][1:] #don't take 0th elmnt; it is the binnum

            interacting_bin_vectors.append(bin_data)
            interacting_bin_freqs.append(freq)

        #if no data vectors to combine, skip it
        assert interacting_bin_vectors != []

        print len(interacting_bin_vectors), gene #display num intxng bins/name

        #to weight intxtr vectors, all iced freqs in effect
        #scaled between 0-1 (by appending them below)
        #goal is to make columns comparable across samples
        # taking average is bad, see eg below
    #    In [34]: q = [[44, 10], [22, 20], [15, 50]] # bin interactors
    #    In [37]: np.average(q, axis=0, weights=[10,20,30])
    #    Out[37]: array([ 22.16666667,  33.33333333])
    #    VS
    #    In [39]: q = [[44, 50], [20, 25], [10, 45]] #similar bin inxts
    #    In [40]: np.average(q, axis=0, weights=[3, 2, 1]) #lower freqs
    #    Out[40]: array([ 30.33333333,  40.83333333  ])
    #    notice low weight vector but corresp high end vals
    #    taking dot fits better with intuition about how
    #    itneracting bins should be combined
    #    VS
    #    In [57]: wqm.dot(q)
    #    Out[57]: array([ 2.76923077,  3.20512821])

        #insert global min and max freq to inxbin vectors
        interacting_bin_freqs.append(min_freq)
        interacting_bin_freqs.append(max_freq)

        interacting_bin_freqs = np.array(interacting_bin_freqs)
        interacting_bin_freqs = interacting_bin_freqs.reshape((interacting_bin_freqs.shape[0],1))
        mms = MinMaxScaler(feature_range=(0,1))
        interacting_bin_freqs = mms.fit_transform(interacting_bin_freqs)
        interacting_bin_freqs = interacting_bin_freqs.ravel()
        interacting_bin_freqs = interacting_bin_freqs[:-2] #remove artificially added values

    #    r = np.average(interacting_bin_vectors, axis=0, weights=interacting_bin_freqs)
        r = interacting_bin_freqs.dot(interacting_bin_vectors)
        sq = sparse.csr_matrix(r) #improve memory footprint by converting to sparse matrix
        print sq.shape, '*'
        assert sq.shape[1] == (precomputed_bins.shape[1] - 1)
        time = datetime.datetime.now() - start
        print time
        log_file = open('./log_file_constructionfinal.txt', 'a')
        log_file.write(gene+'\n')
        log_file.close()
        #return the sparse vector
        
        return (gene, sq)
    except:
        #write error genes ot disk
        log_file = open('./error_log_constructfinal.txt', 'a')
        log_file.write(gene+'\n')
        log_file.close()
        return (gene, None)

def mp_handler():
    p = multiprocessing.Pool(12) #use four processors
    result = p.map(mp_worker, genes)
    joblib.dump(result, './mp_handler_result3')
    return '', ''
#    print result[0][1].shape, '***'
    r = sparse.csr_matrix(np.ones(precomputed_bins.shape[1] - 1)) #seed final matrix with a row of ones; actual size, since binnums cut out

    #once above line complete, vertically stack the vectors
    gene_index = []
    for gene, sp_vector in result:
        if sp_vector == None:
            continue
        print gene, sp_vector.shape
        #before stacking, densify
        row = sp_vector.todense()
        r = sparse.vstack([r, row])
        gene_index.append(gene)

    r = r.tocsc()[1:] #trim seed ones row
    assert r.shape[0] == len(gene_index)
    return r, gene_index


if __name__ == '__main__':
#    genes = get_gene_list('/home/rm/transcription_features/validation_idea/data/RNASeq/EPAS1_siRNA_fpkm/gene_exp.diff')
#    genes = get_gene_list('/home/rm/transcription_features/validation_idea/data/RNASeq/imr90_tnfa/diff_out/gene_exp.diff')
    genes = get_gene_list('./huvec_rnaseq_analysis/DE_genes_ALLGENES.csv')
    print len(genes), 'map to ensemble transcripts'
    #start =datetime.datetime.now()
#    genes = ['c10orf55', 'cc2d1b', 'camk2g', 'alkbh8', 'ttty15']
#    random.shuffle(genes)
    
    mainstart =datetime.datetime.now()
    procstart = mainstart.strftime('%H:%M:%S')
    log_file = open('./timestamp.txt', 'a')
    log_file.write(procstart)
    log_file.close()
    
    r, gene_index = mp_handler()
    joblib.dump(r, './huvec_precomputed/huvec_final_matrix')
    joblib.dump(gene_index, './huvec_precomputed/huvec_final_matrix_genes')
    print r
    print r.shape
