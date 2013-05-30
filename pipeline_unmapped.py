#!/usr/bin/python

import os;
import sys;
import pickle;
import csv;

from Bio import SeqIO;
from Bio import Entrez;

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

def BLAST(query_file, dbloc, outfile, blast_opts=""):

  prog = "blastn";

  opts = "qseqid sseqid slen length mismatch gapopen pident evalue bitscore";
  cmd  = "%s -query '%s' -db '%s' -out '%s' -outfmt '6 %s' %s" % (prog, query, dbloc, outfile, opts, blast_opts);

  return run_cmd(cmd);

#edef

###############################################################################

if C.__unmapped_blast_select_by__ == None or C.__unmapped_blast_select_by__ < 0: 
  print "Please set __unmapped_blast_select_by__ appropriately.";
  sys.exit(1);
#fi

Entrez.email = C.__pipeline_email__;

BO = C.__unmapped_blast_output__();

for i in xrange(len(BO)):
  result = BO[i];
  sn = C.sample_names[i];
  fh = {};
  fo = {};

  print "Finding contaminations for sample %s" % sn;

  print "Parsing BLAST output"; sys.stdout.flush();
  with open(result, 'r') as csvfile:
    r = csv.reader(csvfile, delimiter='\t', quotechar='"');
    for hit in r:
      if (hit[0] in fh and fh[hit[0]][1] > float(hit[C.__unmapped_blast_select_by__])) or (hit[0] not in fh):
        fh[hit[0]] = (hit[1].split('|')[1], float(hit[C.__unmapped_blast_select_by__]));
      #fi
    #efor
  #ewith

  gi_list = [ fh[k][0] for k in fh.keys() ];

  print "Determining source organism of hits"; sys.stdout.flush();
  gi_str = ','.join(gi_list);

  handle  = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb");
  records = SeqIO.parse(handle, "gb");

  k = 0;
  for r in records:
    org = r.annotations['organism'];
    if org in fo:
      fo[org] = (fo[org][0] + 1, fo[org][1] + [ fh[fh.keys()[k]] ] );
    else:
      fo[org] = (1, [ fh[fh.keys()[k]] ]);
    #fi
    k = k + 1;
  #efor

  ll   = [ (k[0], k[1][0]) for k in sorted(fo.items(), key=lambda x: x[1][0], reverse=True) ];
  maxl = max([ len(o[0]) for o in ll ]) + 3;

  print "Possible contaminants found in sample %s:" % sn;
  print "  %-*s %-5s" % (maxl, "Organism", "Hits");
  for o in ll:
    print "  %-*s %-5d" % (maxl, o[0], o[1]);
  #efor
  print "\n\n";
  sys.stdout.flush();

  pickle.dump(fo, open('%s/%s.unmapped_orgs.dat' % (C.outdir, sn), 'w'));
#efor

sys.exit(0);
