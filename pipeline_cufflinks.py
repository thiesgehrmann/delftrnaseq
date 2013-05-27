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

cmd = ("cufflinks -o %s " % C.outdir) + \
      ("%s " % C.cufflinks_opts) + \
      ("-G %s " % C.genome_annot) if C.genome_annot else ("-g %s " % C.genome_guide)  + \
      ("-v ") + \
      ("-b %s " % C.cufflinks_bias_corr if C.cufflinks_bias_corr else "") + \
      ("%s" % C.pre_cufflinks_output());

stdout = open(C.outdir + '/cufflinks_stdout.log', 'w');
stderr = open(C.outdir + '/cufflinks_stderr.log', 'w');

print cmd

retval = run_cmd(cmd, stderr=stderr, stdout=stdout);

stdout.close();
stderr.close();

sys.exit(retval);

