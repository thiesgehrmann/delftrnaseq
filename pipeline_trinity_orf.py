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

ASM = C.__trinity_output__();

prg = "/opt/insy/env/sys_enhance/opt/trinity/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl -t ${file}_trinity.fasta";
pwd = os.getcwd();

for i in xrange(len(ASM)):
  a  = ASM[i];
  sn = C.sample_names[i];

  cmd = prg + \
        " -t " + a;

  run_cmd("mkdir -p %s/trinity_orf" % C.outdir);
  os.chdir("%s/trinity_orf" % C.outdir);

  cmds = [];
  cmds.append(cmd);
  cmds.append("mv 'best_candidates.eclipsed_orfs_removed.cds' '%s/%s.trinity_orfs.fasta'" % (C.outdir, C.sample_names[i]));

  print "Finding ORFs for sample %s." % sn;
  sys.stdout.flush();

  retval = run_seq_cmds(cmds);
  if retval != 0:
    sys.exit(retval);
  #fi

  os.chdir(pwd);
  run_cmd("rm -rf '%s/trinity_orf' % C.outdir");

#efor

sys.exit(0);

