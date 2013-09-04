#!/usr/bin/env python

import os;
import sys;
from ibidas import *
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

L = ",".join([filepair[0] for filepair in TR])
R = ",".join([filepair[1] for filepair in TR])

print "Starting splice site discovery"; sys.stdout.flush();

cmd = "STAR %s --genomeDir %s --genomeLoad LoadAndRemove --readFilesIn %s %s --outSAMstrandField intronMotif --outSAMattributes All --outFilterIntronMotifs RemoveNoncanonical" % (cor(C.star_preal_opts), C.__star_pregg_output__(), L, R);
pid, retval = run_cmd(cmd);
if retval != 0:
   sys.exit(retval);

#splice site filtering
x = Read('SJ.out.tab').Detect()
x = x/('chromosome','intron_start','intron_stop', 'strand', 'motif', 'annotated', 'unique_reads','non_unique_reads', 'max_overhang')
x = x[((_.max_overhang >= C.splice_sites_minmax_overhang) & ((_.unique_reads >= 2) | (_.non_unique_reads >= cor(C.splice_sites_nonunique_reads)))) | (_.unique_reads > cor(C.splice_sites_unique_reads))]
x = x.Get(_.chromosome, _.intron_start, _.intron_stop, _.strand.Each(lambda x: '+' if x==1 else '-',dtype=bytes))
Save(x,'splice_junction_db.tsv')

retval = run_shell("sed 1d splice_junction_db.tsv > splice_junction_db_nohead.tsv")
if retval != 0:
   sys.exit(retval);

cmds = [];
cmds.append("mv Log.out star_prealign.log");
cmds.append("mv Log.progress.out star_prealign.progress.log");
cmds.append("mv Log.final.out star_prealign.final.log");
cmds.append("mv SJ.out.tab star_prealign.SJ.tab");

retval = run_seq_cmds(cmds);
if retval != 0:
   sys.exit(retval);
#fi
  
sys.exit(0);

