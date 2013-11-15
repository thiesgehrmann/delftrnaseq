#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

prg = "Trinity.pl";
UMR = C.__star_al_output_unmapped__();

for i in xrange(len(UMR)):
  l, r = UMR[i];
  sn = C.sample_names[i];

  cmd = prg + \
      (" --seqType fq") + \
      (" --JM 100G") + \
      (" --left %s" % l) + \
      (" --right %s" % r) + \
      (" " + cor(C.trinity_opts));

  cmds = [];
  cmds.append(cmd);
  cmds.append("mv 'trinity_out_dir/Trinity.fasta' '%s/%s.trinity_assembled.fasta'" % (C.outdir, C.sample_names[i]));
  cmds.append("rm -rf trinity_out_dir");

  print "Assembling unmapped reads for sample %s." % sn;
  sys.stdout.flush();

  retval = run_seq_cmds(cmds);
  if retval != 0:
    sys.exit(retval);
  #fi

#efor

sys.exit(0);

