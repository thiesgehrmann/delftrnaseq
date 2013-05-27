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

ASM = C.trinity_output();

prg = "~/env/sys_enhance/opt/trinity/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl -t ${file}_trinity.fasta" 
cmds = [];

for i in xrange(len(ASM)):
  a = ASM[i];

  cmd = prg + \
        " -t " + a;
  cmds.append(cmd);
  cmds.append("mv 'best_candidates.eclipsed_orfs_removed.cds' '%s/%s.trinity_orfs.fasta'" % (C.outdir, C.sample_names[i]));
#efor

print cmds;
#sys.exit(run_seq_cmds(cmds));

