#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

cmds = [];

all_fastq = C.__trimmomatic_output__();

if C.PE:
  S1_sam, S2_sam, S3_sam = C.__isoform_dense_star_align_output_sam__();
  S1_log, S2_log, S3_log = C.__isoform_dense_star_align_output_log__();
  S1_ump, S2_ump, S3_ump = C.__isoform_dense_star_align_output_unmapped__();
else:
  S1_sam = C.__isoform_dense_star_align_output_sam__()[0];
  S1_log = C.__isoform_dense_star_align_output_log__();
  S1_ump = C.__isoform_dense_star_align_output_unmapped__()[0];
#fi

merged_sam           = C.__isoform_dense_star_align_output_merged__();

###############################################################################
  # STAR_ONE

if C.PE:
  star_one_R1 = ','.join([ pair[0] for pair in all_fastq ]);
  star_one_R2 = ','.join([ pair[1] for pair in all_fastq ]);
  RDS = '%s %s' % (star_one_R1, star_one_R2);
else:
  RDS = ','.join(all_fastq);;
#fi

cmd = "STAR %s --genomeDir %s --genomeLoad LoadAndRemove --readFilesIn %s --outSAMattributes All --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped Fastx" % (cor(C.star_al_opts), cor(C.__isoform_dense_genome_generate_dir__), RDS);
cmds.append(cmd);

cmds.append("mv Aligned.out.sam '%s'"  % S1_sam);
cmds.append("mv Log.out '%s'"          % S1_log[0]);
cmds.append("mv Log.progress.out '%s'" % S1_log[1]);
cmds.append("mv Log.final.out '%s'"    % S1_log[2]);
cmds.append("mv SJ.out.tab '%s'"       % S1_log[3]);
if C.PE:
  cmds.append("mv Unmapped.out.mate1 '%s'" % S1_ump[0]);
  cmds.append("mv Unmapped.out.mate2 '%s'" % S1_ump[1]);
else:
  cmds.append("mv Unmapped.out.mate1 '%s'" % S1_ump);
#fi


###############################################################################
  # STAR_TWO
if C.PE:
  cmd = 

  cmd = "STAR %s --genomeDir %s --genomeLoad LoadAndRemove --readFilesIn %s --outSAMattributes All --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped Fastx" % (cor(C.star_al_opts), cor(C.__isoform_dense_genome_generate_dir__), S1_ump[0]);
  cmds.append(cmd);
  cmds.append("mv Aligned.out.sam '%s'"  % S2_sam);
  cmds.append("mv Log.out '%s'"          % S2_log[0]);
  cmds.append("mv Log.progress.out '%s'" % S2_log[1]);
  cmds.append("mv Log.final.out '%s'"    % S2_log[2]);
  cmds.append("mv SJ.out.tab '%s'"       % S2_log[3]);

  cmd = "STAR %s --genomeDir %s --genomeLoad LoadAndRemove --readFilesIn %s --outSAMattributes All --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outReadsUnmapped Fastx" % (cor(C.star_al_opts), cor(C.__isoform_dense_genome_generate_dir__), S1_ump[1]);
  cmds.append(cmd);
  cmds.append("mv Aligned.out.sam '%s'"  % S3_sam);
  cmds.append("mv Log.out '%s'"          % S3_log[0]);
  cmds.append("mv Log.progress.out '%s'" % S3_log[1]);
  cmds.append("mv Log.final.out '%s'"    % S3_log[2]);
  cmds.append("mv SJ.out.tab '%s'"       % S3_log[3]);
#fi



retval = run_seq_cmds(cmds);
if retval != 0:
  sys.exit(retval);
#fi

if C.PE:
  cmd = 'cat %s %s %s > %s' % ( S1_sam, S2_sam, S3_sam, merged_sam);
else:
  cmd = 'cp %s %s' % (S1_sam, merged_sam);
#fi

sys.exit(run_cmd(cmd, shell=True));

