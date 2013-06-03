#!/usr/bin/python

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

bamin  = C.__pre_cufflinks_merge_output__();
bamout = C.__pre_cufflinks_merge_output__();

cmd = 'bamtools sort -byname -in %s -out %s' % (bamin, bamout);

sys.exit(run_cmd(cmd));

