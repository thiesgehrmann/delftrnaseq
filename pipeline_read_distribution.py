#!/usr/bin/env python

from pipeline_common import *
import sys
import HTSeq
import pylab
import matplotlib.pyplot as plot
import numpy as np
import scipy.stats

def read_bam_file(filename, figure_filename, PE, n_bins = 50) :
    reader = HTSeq.BAM_Reader(filename)
    left_read_length, right_read_length, insert_size, read_length = [], [], [], []
    k = 0
    for alg in reader :
        if alg.not_primary_alignment or alg.failed_platform_qc :
            continue
        length = len(alg.read.seq)
        if alg.paired_end :
            k = k + 1
            if alg.proper_pair and alg.mate_aligned and alg.pe_which == 'second' :
                insert_size.append(abs(alg.inferred_insert_size))
            if alg.pe_which == 'first' :
                left_read_length.append(length)
            elif alg.pe_which == 'second' :
                right_read_length.append(length)
        else :
            read_length.append(length)
        if k > 100000 :
            break
    
    # Plot stuff here...
    if PE :
        insert_size = np.array(insert_size)
        values = scipy.stats.mstats.mquantiles(insert_size, prob = [0.01, 0.99])
        values = values.astype(np.int32)
        insert_size = insert_size[insert_size >= values[0]]
        insert_size = insert_size[insert_size <= values[1]]
        f, ax = plot.subplots(1, 3)
        ax[0].set_xlabel('Insert size', fontsize = 16)
        ax[0].set_ylabel('Fraction of pairs', fontsize = 16)
        ax[0].set_xlim(values[0], values[1])
        #n, bins, patches = ax[0].hist(insert_size, n_bins, normed = True, histtype = 'stepfilled')
        n, bins, patches = ax[0].hist(insert_size, bins = range(values[0], values[1] + 1), normed = True, histtype = 'stepfilled')
        pylab.setp(patches, 'facecolor', 'b', 'alpha', 0.45)
        
        ax[1].set_xlabel('Read 1 length', fontsize = 16)
        ax[1].set_ylabel('Fraction of reads', fontsize = 16)
        ax[1].set_xlim(np.min(left_read_length), np.max(left_read_length))
        #n, bins, patches = ax[1].hist(left_read_length, n_bins, normed = True, histtype = 'stepfilled')
        n, bins, patches = ax[1].hist(left_read_length, bins = range(values[0], values[1] + 1), normed = True, histtype = 'stepfilled')
        pylab.setp(patches, 'facecolor', 'b', 'alpha', 0.45)
        
        ax[2].set_xlabel('Read 2 length', fontsize = 16)
        ax[2].set_ylabel('Fraction of reads', fontsize = 16)
        ax[2].set_xlim(np.min(right_read_length), np.max(right_read_length))
        #n, bins, patches = ax[2].hist(right_read_length, n_bins, normed = True, histtype = 'stepfilled')
        n, bins, patches = ax[2].hist(right_read_length, bins = range(values[0], values[1] + 1), normed = True, histtype = 'stepfilled')
        pylab.setp(patches, 'facecolor', 'b', 'alpha', 0.45)
        f.set_size_inches(12, 3)
    else :
        f, ax = plot.subplots(1, 1)
        ax.set_xlabel('Read length', fontsize = 16)
        ax.set_ylabel('Fraction of reads', fontsize = 16)
        ax.set_xlim(np.min(read_length), np.max(read_length))
        #n, bins, patches = ax.hist(read_length, n_bins, normed = True, histtype = 'stepfilled')
        n, bins, patches = ax.hist(read_length, bins = range(value[0], value[1] + 1), normed = True, histtype = 'stepfilled')
        pylab.setp(patches, 'facecolor', 'b', 'alpha', 0.45)
        f.set_size_inches(4, 3)
    f.tight_layout()
    plot.savefig(figure_filename, bbox_inches = 'tight')

C = init_conf()
sort_outputs = C.__post_star_al_sort_output__()
for bam_file in sort_outputs :
    prefix, ext = os.path.splitext(bam_file)
    figure_file = prefix + '.read_stats.pdf'
    print 'Processing file: %s' % bam_file
    read_bam_file(bam_file, figure_file, C.PE)

sys.exit(0)
