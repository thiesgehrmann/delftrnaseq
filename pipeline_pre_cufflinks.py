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

bamout = '%s/bamout.bam' %C.outdir;

cmds = [ 'bamtools merge %s -out %s' % (' -in '.join([''] + C.star_al_output()), bamout),
         'bamtools sort -in %s -out %s' % (bamout, C.pre_cufflinks_output()),
         'rm %s' % (bamout) ];

print cmds

sys.exit(run_seq_cmds(cmds));


