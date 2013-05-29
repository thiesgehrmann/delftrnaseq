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

ALN = C.star_al_output_sam();

cmds_v = [];
cmds_s = []

for i in xrange(len(ALN)):
  a  = ALN[i];
  sn = C.sample_names[i];

  nfname = '%s/%s.star_align.bam' % (C.outdir, sn);
  sfpref = '%s/%s.star_align_sort_name' % (C.outdir, sn);

  cmds_v.append("samtools view -Sb '%s' -o '%s'" % (a, nfname));
  cmds_s.append("samtools sort -n -m 100000000000 '%s' '%s'" % (nfname, sfpref));

#efor

print "Converting SAM to BAM"; sys.stdout.flush();
retval = run_par_cmds(cmds_v, max_threads=C.max_threads);
if retval != 0:
 sys.exit(retval);
#fi 

print "Sorting BAM files by name"; sys.stdout.flush();
retval = run_par_cmds(cmds_s, max_threads=C.max_threads);
if retval != 0:
  sys.exit(retval);
#fi

sys.exit(0);

