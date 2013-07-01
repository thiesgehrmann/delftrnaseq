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

bf  = C.__post_star_al_sort_output__();
gtf = C.__cufflinks_output__()[2];

cmds = [];

for (a,b) in C.cuffdiff_cmp:
  a_reps = [ bf[i] for i in xrange(len(bf)) if C.sample_labels[i] == a ];
  b_reps = [ bf[i] for i in xrange(len(bf)) if C.sample_labels[i] == b ];
  

  cmd = "cuffdiff " \
      + "%s " % cor(C.cuffdiff_opts) \
      + "-L '%s,%s' " % (a, b) \
      + "-b '%s' " % C.genome \
      + "%s " % gtf \
      + "%s " % ','.join(a_reps) \
      + "%s " % ','.join(b_reps);

  cmds.append(cmd);

#efor

#retval = run_seq_cmds(cmds);

print cmds[0];


sys.exit(retval);
