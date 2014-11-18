#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

TR = C.__trimmomatic_output__();



SAMS = cor(C.__star_al_output_sam__);
UNMP = cor(C.__star_al_output_unmapped__);
LOGS = cor(C.__star_al_output_logs__);

for i in xrange(len(TR)):
  print "Starting work on sample '%s'" % C.sample_names[i]; sys.stdout.flush();

  sname = '%s/%s' % (C.outdir, C.sample_names[i]);

  if C.PE:
    RS = ' '.join(TR[i]);
  else:
    RS = TR[i];
  #fi
  cmd = "STAR %s --genomeDir %s --genomeLoad LoadAndRemove --readFilesIn %s --outSAMattributes All --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped Fastx" % (cor(C.star_al_opts), C.__star_gg_output__(), RS);
  
  run_cmd(cmd);

  cmds = [];

  cmds.append("mv Aligned.out.sam '%s'"  % SAMS[i]);
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

