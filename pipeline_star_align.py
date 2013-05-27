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

TR = C.trimmomatic_output();

for i in xrange(len(TR)):
  L, R = TR[i];

  sname = '%s/%s' % (C.outdir, C.sample_names[i])

  cmd = "STAR %s --genomeDir %s --genomeLoad LoadAndRemove --readFilesIn %s %s --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped Fastx" % (C.star_al_opts, C.star_gg_output(), L, R);
  print cmd
  run_cmd(cmd);
  cmd = "samtools view -Sb Aligned.out.sam"
  with open('Aligned.out.bam', 'w') as output_f:
    retval = run_cmd(cmd, stdout=output_f);
    output_f.close();
  #ewith

  if retval != 0:
    sys.exit(retval);
  #fi

  cmds = [];
  cmds.append("samtools sort -m 100000000000 Aligned.out.bam aligned_sort_name");
  cmds.append("samtools sort -n -m 100000000000 Aligned.out.bam aligned_sort_chr");

  cmds.append("rm Aligned.out.sam Aligned.out.bam");
  cmds.append("mv aligned_sort_name.bam '%s.sort_name.bam'" % sname );
  cmds.append("mv aligned_sort_chr.bam  '%s.sort_chr.bam'"  % sname );
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
