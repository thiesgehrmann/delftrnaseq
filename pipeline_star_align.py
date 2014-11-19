#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

TR = C.__trimmomatic_output__();



SAMS = cor(C.__star_al_output_sam__);
BAMS = cor(C.__star_al_output_bam__);
UNMP = cor(C.__star_al_output_unmapped__);
LOGS = cor(C.__star_al_output_logs__);

cmds = [];

for i in xrange(len(TR)):
  print "Starting work on sample '%s'" % C.sample_names[i]; sys.stdout.flush();

  sname = '%s/%s' % (C.outdir, C.sample_names[i]);

  sname = '%s/%s' % (C.outdir, C.sample_names[i])
  tmpdir = '%s/star_tmp' % C.outdir
  if C.PE:
    RS = ' '.join(TR[i]);
  else:
    RS = TR[i];
  #fi

    # Run the BAM version
  cmd = "STAR %s --genomeDir %s --outTmpDir %s --genomeLoad LoadAndRemove --readFilesIn %s --outSAMattributes All --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000" % (cor(C.star_al_opts), C.__star_gg_output__(), tmpdir, RS);
  cmds.append(cmd);

  if C.isoform_dense_analysis::
      # AND the SAM version
    cmd = "STAR %s --genomeDir %s --genomeLoad LoadAndRemove --readFilesIn %s --outSAMattributes All --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped Fastx" % (cor(C.star_al_opts), cor(C.__star_gg_output__), RS); 
    cmds.append(cmd);
    cmds.append("mv Aligned.out.sam '%s'"  % SAMS[i]);
  #fi
  cmds.append("mv Aligned.sortedByCoord.out.bam '%s.star_align.bam'" % sname);
  cmds.append("mv Log.out '%s'"          % LOGS[0][i]);
  cmds.append("mv Log.progress.out '%s'" % LOGS[1][i]);
  cmds.append("mv Log.final.out '%s'"    % LOGS[2][i]);
  cmds.append("mv SJ.out.tab '%s'"       % LOGS[3][i]);

  if C.PE:
    cmds.append("mv Unmapped.out.mate1 '%s'" % UNMP[i][0]);
    cmds.append("mv Unmapped.out.mate2 '%s'" % UNMP[i][1]);
  else:
    cmds.append("mv Unmapped.out.mate1 '%s'" % UNMP[i])
  #fi

  retval = run_seq_cmds(cmds);
  if retval != 0:
    sys.exit(retval);
  #fi
  
#efor

sys.exit(0);

