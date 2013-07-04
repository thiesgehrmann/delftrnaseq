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

BAM  = C.__post_star_al_bam_output__();
SORT = C.__post_star_al_sort_output__();

cmds = []

for i in xrange(len(BAM)):
  b  = BAM[i];
  s  = SORT[i][:-4];

  cmds.append("samtools sort -m 100000000000 '%s' '%s'" % (b, s));

#efor

print "Sorting BAM files"; sys.stdout.flush();
retval = run_par_cmds(cmds, max_threads=C.__max_threads__);
if retval != 0:
  sys.exit(retval);
#fi

sys.exit(0);

