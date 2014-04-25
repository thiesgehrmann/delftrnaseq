#!/usr/bin/env python

from pipeline_common import *
import sys
import numpy as np

C = init_conf()
bam_outputs = C.__post_star_al_bam_output__()
unmapped_outputs = C.__star_al_output_unmapped__()

max_threads = min(C.__max_threads__, 3)
cmd = 'fastqc -o %s --extract -f bam --threads %d -k 10 %s' % (C.outdir, max_threads, ' '.join(bam_outputs))
retval = run_cmd(cmd)
if retval != 0 :
    sys.exit(retval)

unmapped_outputs_str = ' '.join([' '.join(pair) for pair in unmapped_outputs])
cmd = 'fastqc -o %s --extract -f fastq --threads %d -k 10 %s' % (C.outdir, max_threads, unmapped_outputs_str)
retval = run_cmd(cmd)
sys.exit(retval)
