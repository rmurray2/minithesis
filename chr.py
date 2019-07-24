
from sklearn.externals import joblib

#motif_map = open('/home/rm/transcription_features/validation_idea/data/motifmap/HUMAN_hg19_BBLS_1_00_FDR_0_10.bed').read().split('\n')
motif_map = open('HUMAN_hg19_BBLS_1_00_FDR_0_10.bed').read().split('\n')
print len(motif_map)


#start,stop, name,  bbls (1)
chr_d= {}
for i in motif_map[1:-2]:
    

    i = i.split('\t')
    #human:
#        0: chrom
#        1: start
#        2: end
#        3: name
#        4: bbls
#        5: strand
    #mouse:
#        chrom: 14
#        start: 11
#        end: 3
#        name: 7_12
#        bbls: 1
#        strand: 5

    if i == []: continue
    if 'random' not in i[0]:

        if i[0] not in chr_d:
            
            chr_d[i[0]] = [[int(i[1]),int(i[2]), i[3], float(i[4]), i[5]]]
        else:
            chr_d[i[0]].append([int(i[1]),int(i[2]), i[3], float(i[4]), i[5]])

for chrom, sites in chr_d.iteritems():
    sites.sort()

#joblib.dump(chr_d, '/home/rm/transcription_features/validation_idea/data/motifmap/sites/motif_map_sites_hg19')
joblib.dump(chr_d, 'motif_map_sites_hg19')

