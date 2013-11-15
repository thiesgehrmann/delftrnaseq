#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

cmds = [];

for i in xrange(len(C.samples)):
  L, R = C.samples[i]
  cmds = cmds + [C.trimmomatic_prg + ' %s -phred33 -trimlog %s/%s_trimmomatic.log %s %s %s %s %s %s %s' % 
                       (cor(C.trimmomatic_opts), C.outdir, C.sample_names[i], L, R, 
                       C.outdir + '/' + C.sample_names[i] + '_R1.fastq',  '/dev/null',
                       C.outdir + '/' + C.sample_names[i] + '_R2.fastq', '/dev/null', C.trimmomatic_trim) ];


sys.exit(run_seq_cmds(cmds));

