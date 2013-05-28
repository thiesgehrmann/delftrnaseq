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

for i in xrange(len(ALN)):
  a  = ALN[i];
  sn = C.sample_names[i];

  nfname = '%s/%s.star_align.bam' % (C.outdir, sn);
  sfpref = '%s/%s.star_align_sort_name' % (C.outdir, sn);

  cmds = [];

  cmds.append("samtools view -Sb '%s' -o '%s'" % (a, nfname));
  cmds.append("samtools sort -n -m 100000000000 '%s' '%s'" % (nfname, sfpref));

  print "Sorting sample %s" % sn;
  sys.stdout.flush();

  retval = run_seq_cmds(cmds);
  if retval != 0:
    sys.exit(retval);
  #fi

#efor

retval = run_seq_cmds(cmds);
if retval != 0:
  sys.exit(retval);
#fi

sys.exit(0);

