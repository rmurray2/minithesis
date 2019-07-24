from __future__ import division
import numpy as np
import itertools
import datetime
import multiprocessing
from scipy import sparse
import overlap #built with cython to optimize speed; makes calls to this func lightning fast
from sklearn.externals import joblib
import pandas as pd

################
#load parsed motifmap file
print 'loading parsed motifmap file'
chr_q = joblib.load('./motif_map_sites_hg19') #load parsed flat prediction file
#chr_q = {'chr'+str(k) : v for k, v in chr_q.iteritems()}
################

#################
## build unique TF list
print 'building unique tf lsit'
all_tfs = []

df = pd.read_csv('HUMAN_hg19_BBLS_1_00_FDR_0_10.bed', delimiter='\t', names=range(0, 6))

for index, row in df.iterrows():
    all_tfs.append(row[3])

all_tfs = list(set(all_tfs))
all_tfs.sort()
print (all_tfs)
print len(all_tfs)

#load parsed chipseq data
print 'loading chipseq data'
#chipseq_data = joblib.load('./processed_chipseq/chipseq_dict2')
chipseq_data = joblib.load('./processed_chipseq/chipseq_dict')

#################

#ideally this would be a script to generate sparse mat of precomputed individual bin vectors from inputs:
    #main chipseq dir containing peak files
        #prog would call preprocess_chipseq.py
    #bed file containing tfbs predictions
        #prog would parse file
    #HiC-pro _abs.bed file
    #output sparse matrix fn; output feature labels list fn

def find_sites(ch, start, stop, ):
    '''
    finds sites in motifmap http://motifmap.ics.uci.edu hg19 precomputed sites
    '''

    if 'chr' not in ch:
        ch = 'chr'+str(ch)
    doi = chr_q[ch] #zoom in in chr of interest
    result = []
    hit_patch = False
    null = False
    for i in doi: #for each TF site on that chr:
        site_start, site_stop, name, bbls = i[0], i[1], i[2], i[3]

        if overlap.is_overlapping(start, stop, site_start, site_stop):
            result.append([site_start, site_stop, name, bbls])
            hit_patch = True
            null = True
        elif hit_patch == True: # if it doesn't overlap and a patch of overlap has already been hit,
                                #then break out instead of looping unnecessarily
            break

    return result

def calc_distance(tf_motif, chrom_mark):
    #tf_i_motifs is of form [start, end, bbls]
    #chrom_mark is of form [start, end, signalvalue, summit]
    #produces a distance metric which is always between 0-1
    #values close to 1 mean the tfbs is close to the chipseq summit
    #lower values mean the tfbs is toward one of the peak boundaries

    summit = chrom_mark[-1]
    tfbs_midpoint = int((abs(tf_motif[0] - tf_motif[1]) / 2) + tf_motif[0])

    numerator = abs(summit - tfbs_midpoint) #distance between chrom mark summit and tbfs
    denominator = abs(chrom_mark[0] - chrom_mark[1]) / 2 #half length of chromatin mark
    fraction = numerator / denominator
    return 1 - fraction

#pseudocode
#dv for each bin is of form:
#    [#TF occurrances in bin
#     max bbls of TFs
#        #TFs in mark
#        max signal value among marks containing TF
#        max (distance measure) (1 - ( |summit - TFBSloc| / len( 1/2 mark) ) )
#
#((5 marks * 3 data points/mark) * 607 TFs) + 607 TFs X 2
#+
#whatever the self data vector ends up being
#
#d = {}
#for each bin:
#    get sites given bin range
#    for each chipseq mark
#        get all marks overlapping with bin
#        {mark : o[site...]]
#      = []
#    for each tf
#        inner_temp= []
#        append #TF occurrances to inner_temp
#        get max bbls of all TFs in bin
#        for each mark in mark list:
#            get #TFs which overlap with marks
#            get max (distance measure) (1 - ( |summit - TFBSloc| / len( 1/2 mark) ) )
#            get max signal value among marks containing tf
#
#            append 4 data points to inner temp
#        add innertemp to temp
#    d[bin#] = temp

def precompute_bin_vectors(binfile):
#    r = sparse.csr_matrix(np.ones(10320)) #len of vect plus one
    n_chipseq_datapoints = 3
    n_tf_specific_datapoints = 2
    n_data_cols = ((len(chipseq_data.keys()) * n_chipseq_datapoints) + n_tf_specific_datapoints) * len(all_tfs)
    n_data_cols = n_data_cols + 1 #add one to len

    r = sparse.csr_matrix(np.ones(n_data_cols))
    count = 0
    with open(binfile) as infile:

        count = 0
        for line in infile:
            count += 1
#            if count == 500:
#                break
#            print count
            if 'chrM' in line or line == '':
                continue
            quit = False
            sline = line.split('\t')

            #find all motifmap sites in specified range (of the bin)
            chrom, start, end, binnum = sline[0], int(sline[1]), int(sline[2]), int(sline[3])
            sites = find_sites(chrom, int(start), int(end))

            #make a motif dict for the bin
            #{tf name : [start, stop, bbls]}
            bin_motifs = {}
            for s in sites:
                if s[2] not in bin_motifs:
                    bin_motifs[s[2]] = [[s[0], s[1], s[3]]]
                else:
                    bin_motifs[s[2]].append([s[0], s[1], s[3]])

            #write to log file
            log_file = open('./log_file.txt', 'a')
            log_file.write(line+'\n')
            log_file.close()
            if sites == []:
                continue

            #find all chipseq marks in bin
            # {experimentname : [ [start, end, signalvalue, summit] ]  }
            marks = {i : [] for i in chipseq_exps}
            for exp in chipseq_exps: #for experiment in sorted chipseq experiment list
                if chrom not in chipseq_data[exp]:
                    continue
                coi = chipseq_data[exp][chrom] #chromosome of interest in experiment of itnerest
                for chipseq_site in coi:
                    if overlap.is_overlapping(chipseq_site[0], chipseq_site[1],start, end): #compare bin coordinates with chipseq mark coordintaes, check for overlap
                        marks[exp].append(chipseq_site)
                        #quit = True #quit program after encountering first bin containing some mark

            #temp is the container for the bin's main data vector
            temp = []
            for tf in all_tfs:   #for each tf in complete tf list

                if tf not in bin_motifs:
                    for i in range ((len(chipseq_exps) * n_chipseq_datapoints) + n_tf_specific_datapoints): #if ith tf motif not present,
                        temp.append(0)                                                   #empty vector
                    continue

                #if control reaches this point,
                #it means the tfbs is present in the bin
                inner_temp = []   #inner temp holds temporary tf-specifc/chipseq-specific data

                inner_temp.append(len(bin_motifs[tf]))  #total number of tf motifs in bin
                inner_temp.append(max([j[-1] for j in bin_motifs[tf]]))  #max bbls among tf motifs in bin

                for exp in chipseq_exps:  #for each chip experiment
                    #the  purpose of this loop is to fill up chipseq_temp
                    #by the end of each iteration, chipseq_temp should have 3 data points
                    #(3 data points for each chip seq)
                    chipseq_temp = []
                    if marks[exp] == []: #if chromatin mark 'exp' NOT present in
                                         # in the bin, append empty vect
                        for i in range(n_chipseq_datapoints):
                            chipseq_temp.append(0)

                    else: #if chromatin mark 'exp' IS present in the bin

                        tf_i_motifs_in_mark_j = []
                        mark_j_containing_tf_i = []
                        tfbs_summit_distances = []

                        #for each tfbs, check whether it falls in a chipseq mark
                        for tf_i_motifs in bin_motifs[tf]: #for each ith tfbs in bin
                            #tf_i_motifs is of form [start, end, bbls]
                            tfbs_start, tfbs_end = tf_i_motifs[0], tf_i_motifs[1]
                          #  print tfbs_start, tfbs_end
                            for chrom_mark in marks[exp]: #for each mark from chipseq experiment j
                                #chrom_mark is of form [start, end, signalvalue, summit]
                                mark_start, mark_end = chrom_mark[0], chrom_mark[1]
                                if overlap.is_overlapping(tfbs_start, tfbs_end, mark_start, mark_end):
                                    #print tf_i_motifs, chrom_mark, "***"
                                    tf_i_motifs_in_mark_j.append(tf_i_motifs)
                                    mark_j_containing_tf_i.append(chrom_mark)

                                    #for each overlap if tfbs/chrom mark, calculate distance measure
                                    dm = calc_distance(tf_i_motifs, chrom_mark)
                                    #print 'dm', dm
                                    tfbs_summit_distances.append(dm)

                        #remove dupes
                        mark_j_containing_tf_i.sort()
                        tf_i_motifs_in_mark_j.sort()

                        mark_j_containing_tf_i = list(mark_j_containing_tf_i for mark_j_containing_tf_i,_
                                                      in itertools.groupby(mark_j_containing_tf_i))
                        tf_i_motifs_in_mark_j = list(tf_i_motifs_in_mark_j for tf_i_motifs_in_mark_j,_
                                                      in itertools.groupby(tf_i_motifs_in_mark_j))

                        #tf_i_motifs_in_mark_j is of form [start, end, bbls]
                        #mark_j_containing_tf_i is of form [start, end, signalvalue, summit]

                        #if at least one tfbs overlaps with a chipseq mark
                        #tfbs/chipseq data points can be generated
                        #loop below is the heart of where the tfbs data is related to the chipseq data
                        if len(tf_i_motifs_in_mark_j) > 0:
                            chipseq_temp.append(len(tf_i_motifs_in_mark_j)) # num of ith tf sites which overlap
                                                                            #with mark j sites
                            chipseq_temp.append(max([j[2] for j in mark_j_containing_tf_i]))#get chipseq max signal
                                                                                            #value among
                                                                                            #marks containing tf
                            chipseq_temp.append(max(tfbs_summit_distances)) #among tfbs overlapping with chip
                                                                            #site, the dist of the closest one
                                                        #to peak summit

                        else: #if no tfs in any marks, then append empty vector
                            for i in range(n_chipseq_datapoints):
                                chipseq_temp.append(0)

                    #append tfbs-chipseq data (3 data points)
                    for i in chipseq_temp:
                        inner_temp.append(i)

                #append tfbs-chipseq data to overall vector
                for i in inner_temp:
                    temp.append(i)

            #insert the binnumber as the first element in the bins data vector
            temp.insert(0, binnum)
            sq = sparse.csr_matrix(temp) #improve memory footprint by converting to sparse matrix
            row = sq.todense()
            r = sparse.vstack([r, row])
    r = r.tocsc()[1:]
    return r


if __name__ == '__main__':

    chipseq_exps = chipseq_data.keys()
    chipseq_exps.sort()
    print chipseq_exps
    feature_labels = []
    for tfbs in all_tfs:
        feature_labels.append('num_' + tfbs)
        feature_labels.append('max_bbls_' + tfbs)
        for exp in chipseq_exps:
            feature_labels.append('#_' + tfbs + '_sites_overlapping_'+exp)
            feature_labels.append('max_signal_among_' + exp + '_containing_' + tfbs)
            feature_labels.append('max_' + tfbs + 'summit_distance_' + exp)

    joblib.dump(feature_labels, './motifmap_cutoff_BroadChipSeq_protein_coding_data/huvec_feature_labels')
    print len(feature_labels), 'flabels'

#    bin_file = './HiC/huvec_validation/huvec_results/hic_results/matrix/huvec/raw/30000/huvec_30000_abs.bed'
    bin_file = './huvec_out/hic_results/matrix/huvec/raw/30000/huvec_30000_abs.bed'

    start_time = datetime.datetime.now()
    res = precompute_bin_vectors(bin_file)
    print datetime.datetime.now() - start_time
#    joblib.dump(res, './huvec_precomputed/huvec_precomputed_bins')
    joblib.dump(res, './huvec_precomputed/huvec_precomputed_bins2')
