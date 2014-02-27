#!/usr/bin/env python
from pipeline_common import *;
from ibidas import *

C = init_conf()

TR = C.__trimmomatic_output__();

L = ",".join([filepair[0] for filepair in TR])
R = ",".join([filepair[1] for filepair in TR])

print "Starting splice site discovery"; sys.stdout.flush();

cmd = "STAR %s --genomeDir %s --genomeLoad LoadAndRemove --readFilesIn %s %s --outSAMstrandField intronMotif --outSAMattributes All --outFilterIntronMotifs RemoveNoncanonical" % (cor(C.star_preal_opts), C.__star_pre_splice_output__(), L, R);
retval = run_cmd(cmd);
if retval != 0:
   sys.exit(retval);

#splice site filtering
x = Read('SJ.out.tab').Detect()
x = x/('chromosome','intron_start','intron_stop', 'strand', 'motif', 'annotated', 'unique_reads','non_unique_reads', 'max_overhang')
x = x[((_.max_overhang >= C.splice_sites_minmax_overhang) & ((_.unique_reads >= 2) | (_.non_unique_reads >= cor(C.splice_sites_nonunique_reads)))) | (_.unique_reads > cor(C.splice_sites_unique_reads))]
x = x.Get(_.chromosome, _.intron_start, _.intron_stop, _.strand.Each(lambda x: '+' if x==1 else '-',dtype=bytes))
Save(x,'%s/splice_junction_db.tsv' % C.outdir)

retval = run_shell("sed 1d %s/splice_junction_db.tsv > %s/splice_junction_db_nohead.tsv" % (C.outdir, C.outdir))
if retval != 0:
   sys.exit(retval);

cmds = [];
cmds.append("mv Log.out %s/star_prealign.log" % C.outdir);
cmds.append("mv Log.progress.out %s/star_prealign.progress.log" % C.outdir);
cmds.append("mv Log.final.out %s/star_prealign.final.log" % C.outdir);
cmds.append("mv SJ.out.tab %s/star_prealign.SJ.tab" % C.outdir);

retval = run_seq_cmds(cmds);
if retval != 0:
   sys.exit(retval);
#fi
  
sys.exit(0);

