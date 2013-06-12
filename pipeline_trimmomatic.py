#!/usr/bin/python

import os;
import sys;

from pipeline_common import *;

###############################################################################

def usage(a1):
  print "Usage:  %s <config file>" % a1;
#edef

if len(os.sys.argv) != 2:
  usage(os.sys.argv[0]);
  os.sys.exit(1);
#fi

C = PIPELINECONF(os.sys.argv[1]);
run_cmd('mkdir -p %s' % C.outdir);

###############################################################################

prg="trimmomaticPE"

cmds = [];

for i in xrange(len(C.samples)):
  L, R = C.samples[i]
  cmds = cmds + ['trimmomaticPE %s -phred33 -trimlog %s/%s_trimmomatic.log %s %s %s %s %s %s %s' % (cor(C.trimmomatic_opts), C.outdir, C.sample_names[i], L, R, 
                                                                                           C.outdir + '/' + C.sample_names[i] + '_R1.fastq',  '/dev/null',
                                                                                           C.outdir + '/' + C.sample_names[i] + '_R2.fastq',  '/dev/null', C.trimmomatic_trim) ];
#efor

sys.exit(run_seq_cmds(cmds));

