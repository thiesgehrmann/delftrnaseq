#!/usr/bin/env python
from pipeline_common import *;

import pickle;
import csv;
import time
from ibidas import *

from Bio import SeqIO;
from Bio import Entrez;


C = init_conf()

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

  gi_list = [ fh[k][0] for k in fh.keys() if fh[k][1] > 1500];

  print "Determining source organism of %d hits" % len(gi_list); sys.stdout.flush();

  k = 0;
  xsize = 50;
  for j in xrange(0,len(gi_list), xsize):
    gi_str = ','.join(gi_list[j:(j+xsize)]);

    print "Sending query %d" % (k/xsize); sys.stdout.flush();
    print gi_str;
    r = 0
    handle = None
    while 1:
        try:
            handle  = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb");
        except Exception, e:
            r = r + 1
            time.sleep(60)
            if r> 30: 
                raise
            else:
                print e
                continue
        break


    print "Parsing query"; sys.stdout.flush();

    l = '\n';
    while l != '':
      l = handle.readline();
      if l[2:10] != 'ORGANISM':
        continue;
      #fi
      org = l[10:];
      if org in fo:
        fo[org] = (fo[org][0] + 1, fo[org][1] + [ fh[fh.keys()[k]] ] );
      else:
        print org; sys.stdout.flush();
        fo[org] = (1, [ fh[fh.keys()[k]] ]);
      #fi
      k = k + 1;
    #efor
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
  Save(fo, '%s/%s.unmapped_orgs.dat' % (C.outdir, sn));
#efor

sys.exit(0);
