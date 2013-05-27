#!/usr/bin/python

import os;
import sys;
import pickle;

from Bio.Blast import NCBIWWW;
from Bio.Blast import NCBIXML;
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

#C = PIPELINECONF(os.sys.argv[1]);
#run_cmd('mkdir -p %s' % C.outdir);

###############################################################################

Entrez.email = 'thiesgehrmann@gmail.com' #C.pipeline_email;

#ORF = C.trinity_orf_output();
ORF = [ 'test.fasta' ];
oo = [ {} for x in ORF ];

for i in xrange(len(ORF)):
  fa = ORF[i];
  sn = ''#C.sample_names[i];
  fo = {};

  print "Reading Fasta File"  
  fd = open(fa, 'r');
  S = fd.read();
  fd.close();

  print "Performing BLAST query"
  handle = NCBIWWW.qblast('blastn', 'nr', S);

  print "Parsing BLAST output";
  blastres = NCBIXML.parse(handle);
  gids = [];
  for q in blastres:
    if len(q.alignments) < 1:
      continue;
    #fi
    id = q.alignments[0].hit_id.split('|');
    if len(id) > 1:
      gids.append(q.alignments[0].hit_id.split('|')[1]);
    #fi
  #efor

  print "Determining source organism of hits"
  gi_str = ','.join(gids);
  handle  = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb");
  records = SeqIO.parse(handle, "gb");

  for r in records:
    org = r.annotations['organism'];
    if org in fo:
      fo[org] = fo[org] + 1;
    else:
      fo[org] = 1;
    #fi
  #efor

  print i
  oo[i] = fo;
  print oo[i];

  pickle.dump(oo, open('otherorgs.py', 'w'));
#efor

