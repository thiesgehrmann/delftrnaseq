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

ALN = C.__star_al_output_sam__();
BAM = C.__post_star_al_bam_output__();

cmds = [];

for i in xrange(len(ALN)):
  a = ALN[i];
  b = BAM[i];

  cmds.append("samtools view -Sb '%s' -o '%s'" % (a, b));

#efor

print "Converting SAM to BAM"; sys.stdout.flush();
retval = run_par_cmds(cmds, max_threads=C.__max_threads__);
if retval != 0:
 sys.exit(retval);
#fi 

sys.exit(0);

