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

gffout = C.__genome_annot_format_output__();

cmd = 'gffread ' \
    + '-E %s ' % C.genome_annot \
    + '-o %s'  % gffout;

retval = run_cmd(cmd);

if retval != 0:
    print "ERROR: gffread fails, simply copying genome annotation file instead!!!!!"
    cmd = 'cp %s %s' % (C.genome_annot, gffout)
    retval = run_cmd(cmd)

sys.exit(retval);
