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

def BLAST_cmd(query_file, dbloc, outfile, fields="qseqid sseqid slen length mismatch gapopen pident evalue bitscore", blast_opts=""):

  prog = "blastn";

  cmd  = "%s -query '%s' -db '%s' -out '%s' -outfmt '6 %s' %s" % (prog, query_file, dbloc, outfile, fields, blast_opts);

  return cmd;

#edef

###############################################################################

if C.blast_db == None:
  print "You need to provide a BLAST database. Please set the variable 'blast_db'.";
  sys.exit(1);
#fi
if C.__unmapped_blast_fields__ == None:
  print "You need to specify which fields you want BLAST to return.\nPlease set __unmapped_blast_fields__ appropriately.";
  sys.exit(1);
#fi;
  

###############################################################################

ORF = C.__trinity_orf_output__();

blast_db = ''.join(C.blast_db.split('.')[0:-1]);

cmds = [];

for i in xrange(len(ORF)):
  fa = ORF[i];
  sn = C.sample_names[i];
  outf = '%s/%s.unmapped_blast.tsv' % (C.outdir, sn);

  cmds.append("echo 'Running BLAST for sample %s'." % sn);
  cmds.append(BLAST_cmd(fa, blast_db, outf, C.__unmapped_blast_fields__, cor(C.unmapped_blast_opts)));
#efor


retval = run_seq_cmds(cmds);
if retval != 0:
  sys.exit(retval);
#fi

sys.exit(0);
