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

prg="STAR"

cmds = [ "mkdir -p %s " % C.__star_gg_output__(),
         "%s --runMode genomeGenerate --genomeDir %s %s --genomeFastaFiles %s" % (prg, C.__star_gg_output__(), cor(C.star_gg_opts), ' '.join(C.genome)) ]

sys.exit(run_seq_cmds(cmds));


