from __future__ import division
from sys import argv
from sklearn.externals import joblib
import os
import pandas as pd

#input dir is given, and dict is created
# chr : [ [start, end, signalvalue, summit] ]
# with the value lists sorted

def main(di, outdir):
    d = os.listdir(di)

    od = {}
    #overdict (od) is { chipseq_expname : { chr : [[...]] } }
    for f in d:
        #process narrowpeak and broadpeak files differently 
        #this doesn't look necessary; look getting rid of conditionals
        o = {}
        if 'narrow' in f:
            df = pd.read_table(di + f, header=None)
            for index, row in df.iterrows():
                if row[0] not in o:
                    summit = min(row[1], row[2]) + (abs(row[1] - row[2]) / 2)
                    o[row[0]] = [[row[1], row[2], row[6], int(summit)]]
                else:
                    summit = min(row[1], row[2]) + (abs(row[1] - row[2]) / 2)
                    o[row[0]].append([row[1], row[2], row[6], int(summit)])

            for chrom, sites in o.iteritems():
                sites.sort()
            od[f] = o #
            print (f, 'processed')

        elif 'broad' in f:

            df = pd.read_table(di + f, header=None)
            for index, row in df.iterrows():
                if row[0] not in o:
                    summit = min(row[1], row[2]) + (abs(row[1] - row[2]) / 2)
                    o[row[0]] = [[row[1], row[2], row[6], int(summit)]]
                else:
                    summit = min(row[1], row[2]) + (abs(row[1] - row[2]) / 2)
                    o[row[0]].append([row[1], row[2], row[6], int(summit)])
            for chrom, sites in o.iteritems():
                sites.sort()
            od[f] = o #
            print (f, 'processed')
        else:

            print (f, 'skipped')

    print ('done')
    if outdir[-1] != '/':
        joblib.dump(od, outdir+'/chipseq_dict2')
    else:
        joblib.dump(od, outdir+'chipseq_dict2')


if __name__ == '__main__':
    _, di, outdir = argv
    '''
    di is directory to process containing narrowpeak or broadpeak files
    outdir is the desired output directory of processed chipseq data
    '''
    main(di, outdir)
