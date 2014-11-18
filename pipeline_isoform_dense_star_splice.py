#!/usr/bin/env python
from pipeline_common import *;
from ibidas import *

C = init_conf()

TR = C.__trimmomatic_output__();

if C.PE:
  L = ",".join([filepair[0] for filepair in TR])
  R = ",".join([filepair[1] for filepair in TR])
else:
  L = ','.join(TR);
  R = "";
#fi

LOGS = cor(C.__isoform_dense_star_splice_output_logs__);

print "Starting splice site discovery"; sys.stdout.flush();

cmd = "STAR %s --genomeDir %s --genomeLoad LoadAndRemove --readFilesIn %s %s --outSAMstrandField intronMotif --outSAMattributes All --outFilterIntronMotifs RemoveNoncanonical" % (cor(C.star_preal_opts), cor(C.__isoform_dense_star_pre_splice_output__), L, R);
retval = run_cmd(cmd);
if retval != 0:
   sys.exit(retval);

#splice site filtering
x = Read('SJ.out.tab').Detect()
x = x/('chromosome','intron_start','intron_stop', 'strand', 'motif', 'annotated', 'unique_reads','non_unique_reads', 'max_overhang')
x = x[((_.max_overhang >= C.splice_sites_minmax_overhang) & ((_.unique_reads >= 2) | (_.non_unique_reads >= cor(C.splice_sites_nonunique_reads)))) | (_.unique_reads > cor(C.splice_sites_unique_reads))]
x = x.Get(_.chromosome, _.intron_start, _.intron_stop, _.strand.Each(lambda x: '+' if x==1 else '-',dtype=bytes))
Save(x, cor(C.__isoform_dense_star_splice_output__), names=False);

cmds = [];
cmds.append("mv Log.out %s" % LOGS[0]);
cmds.append("mv Log.progress.out %s" % LOGS[1]);
cmds.append("mv Log.final.out %s" % LOGS[2]);
cmds.append("rm SJ.out.tab");
cmds.append("rm Aligned.out.sam");
cmds.append("rm Unmapped.out.mate1");
if C.PE:
  cmds.append("rm Unmapped.out.mate2");
#fi


retval = run_seq_cmds(cmds);
if retval != 0:
   sys.exit(retval);
#fi
  
sys.exit(0);

