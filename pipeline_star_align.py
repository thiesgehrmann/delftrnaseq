#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

TR = C.__trimmomatic_output__();

for i in xrange(len(TR)):
  print "Starting work on sample '%s'" % C.sample_names[i]; sys.stdout.flush();

  sname = '%s/%s' % (C.outdir, C.sample_names[i])
  tmpdir = '%s/star_tmp' % C.outdir
  if C.PE:
    RS = ' '.join(TR[i]);
  else:
    RS = TR[i];
  #fi
  cmd = "STAR %s --genomeDir %s --outTmpDir %s --genomeLoad LoadAndRemove --readFilesIn %s --outSAMattributes All --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000" % (cor(C.star_al_opts), C.__star_gg_output__(), tmpdir, RS);
  
  run_cmd(cmd);

  cmds = [];

  cmds.append("mv Aligned.sortedByCoord.out.bam '%s.star_align.bam'" % sname);
  cmds.append("mv Log.out '%s.star_align.log'" %sname );
  cmds.append("mv Log.progress.out '%s.star_align.progress.log'" % sname );
  cmds.append("mv Log.final.out '%s.star_align.final.log'" % sname );
  cmds.append("mv SJ.out.tab '%s.star_align.SJ.tab'" % sname );
  cmds.append("mv Unmapped.out.mate1 '%s.star_align_unmapped_R1.fastq'" % sname);
  if C.PE:
    cmds.append("mv Unmapped.out.mate2 '%s.star_align_unmapped_R2.fastq'" % sname);
  #fi

  retval = run_seq_cmds(cmds);
  if retval != 0:
    sys.exit(retval);
  #fi
  
#efor

sys.exit(0);

