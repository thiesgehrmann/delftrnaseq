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

bf  = C.__star_al_output_name__();

gtf = "%s/transcripts.gtf" % C.outdir;


cmd = "cuffdiff " \
    + "-L '%s' " % ','.join([str(l) for l in C.sample_labels]), 
    + "%s " % gtf \
    + "%s " % ' '.join(bf);

sys.exit(run_cmd(cmd));

