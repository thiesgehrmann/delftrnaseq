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

prg = "~/env/sys_enhance/opt/trinity/Trinity.pl";
UMR = C.star_al_output_unmapped();

cmds = [];

for i in xrange(len(UMR)):
  l, r = UMR[i];

  cmd = prg + \
      (" --seqType fq") + \
      (" --JM 100G") + \
      (" --left %s" % l) + \
      (" --right %s" % r) + \
      (" " + C.trinity_opts);

  cmds.append(cmd);

  cmds.append("mv 'trinity_out_dir/Trinity.fasta' '%s/%s_trinity_assembled.fasta'" % (C.outdir, C.sample_names[i]));
#efor

print cmds

#sys.exit(run_seq_cmd(cmds));

