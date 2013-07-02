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

TR = C.__trimmomatic_output__();

for i in xrange(len(TR)):
  L, R = TR[i];

  print "Starting work on sample '%s'" % C.sample_names[i]; sys.stdout.flush();

  sname = '%s/%s' % (C.outdir, C.sample_names[i])

  cmd = "STAR %s --genomeDir %s --genomeLoad LoadAndRemove --readFilesIn %s %s --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped Fastx" % (cor(C.star_al_opts), C.__star_gg_output__(), L, R);

  run_cmd(cmd);

  cmds = [];

  cmds.append("mv Aligned.out.sam '%s.star_align.sam'" % sname);
  cmds.append("mv Log.out '%s.star_align.log'" %sname );
  cmds.append("mv Log.progress.out '%s.star_align.progress.log'" % sname );
  cmds.append("mv Log.final.out '%s.star_align.final.log'" % sname );
  cmds.append("mv SJ.out.tab '%s.star_align.SJ.tab'" % sname );
  cmds.append("mv Unmapped.out.mate1 '%s.star_align_unmapped_R1.fastq'" % sname);
  cmds.append("mv Unmapped.out.mate2 '%s.star_align_unmapped_R2.fastq'" % sname);

  retval = run_seq_cmds(cmds);
  if retval != 0:
    sys.exit(retval);
  #fi
  
#efor

sys.exit(0);

