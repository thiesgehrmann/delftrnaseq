#!/usr/bin/env python

import os;
import sys;
import pickle;
import csv;
import time
from ibidas import *

from Bio import SeqIO
from Bio import Entrez
import MySQLdb as mdb

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

if C.__unmapped_blast_select_by__ == None or C.__unmapped_blast_select_by__ < 0: 
  print "Please set __unmapped_blast_select_by__ appropriately.";
  sys.exit(1);
#fi

Entrez.email = C.__pipeline_email__;

BO = C.__unmapped_blast_output__();
criterion = C.__unmapped_blast_select_by__
cutoff = C.__unmapped_blast_select_by_cutoff__

for i in xrange(len(BO)):
  result = BO[i];
  sn = C.sample_names[i];

  best_hit = {}
  gi2org = {}

  print "Finding contaminations for sample %s" % sn;

  print "Parsing BLAST output"; sys.stdout.flush();
  with open(result, 'r') as csvfile:
    r = csv.reader(csvfile, delimiter='\t', quotechar='"');
    for hit in r:
      query_id = hit[0]
      value = float(hit[criterion])
      gi = hit[1].split('|')[1]
      if not best_hit.has_key(query_id) :
          best_hit[query_id] = ([gi], value)
      elif best_hit[query_id][1] > value :
        best_hit[query_id] = ([gi], value)
      elif best_hit[query_id][1] == value :
        best_hit[query_id][0].append(gi)

  gi_list = []
  for key in best_hit.keys() :
      ids, value = best_hit[key]
      ids = list(set(ids))
      best_hit[key] = (ids, value)
      if value >= cutoff :
          for gi in ids :
              gi_list.append(gi)
  
  gi_list = list(set(gi_list))

  print "[i] Determining source organism of %d hits" % len(gi_list); sys.stdout.flush();

  if not C.unmapped_use_mysql :
    k = 0
    xsize = 50
    for j in xrange(0,len(gi_list), xsize) :
        gi_str = ','.join(gi_list[j:(j+xsize)])

        print "   [i] Sending query %d" % (k/xsize); sys.stdout.flush();
        r = 0
        handle = None
        while True :
            try :
                handle  = Entrez.efetch(db = "nuccore", id = gi_str, rettype = "gb");
            except Exception, e :
                r = r + 1
                time.sleep(60)
                if r > 30 :
                    raise
                else:
                    print e
                    continue
            break


        print "      [i] Parsing query"; sys.stdout.flush();

        l = '\n';
        while l != '':
            l = handle.readline();
            if l[2:10] != 'ORGANISM':
                continue;
            org = l[10:]
            gi = gi_list[k]
            gi2org[gi] = org
            print '         [+] %s -> %s' % (gi, org)
            k = k + 1;

  else :
        n_hits = len(gi_list)
        gi_str = ','.join(gi_list)
        print "   [i] Sending query (size %d)" % n_hits
        con, data = None, None
        try :
            con = mdb.connect(C.unmapped_mysql_host, C.unmapped_mysql_user, C.unmapped_mysql_pass, C.unmapped_mysql_db)
            cur = con.cursor()
            str = 'SELECT gi, source FROM gi2source WHERE gi in (%s)' % gi_str
            cur.execute(str)
            data = cur.fetchall()
        except mdb.Error, e :
            print "Error %d: %s" % (e.args[0], e.args[1])
            sys.stdout.flush()
        finally :
            if con :
                con.close()

        print "      [i] Parsing query results"
        sys.stdout.flush()
        for gi, org in data :
            gi2org[gi] = org
            print '         [+] %s -> %s' % (gi, org)

  cont = {}
  for key in best_hit.keys() :
      gis, score = best_hit[key]
      if score >= cutoff :
          for gi in gis :
              org = gi2org[gi]
              if cont.has_key(org) :
                  hits, ids = cont[org]
                  cont[org] = (hits + 1, ids + [gi])
              else :
                  cont[org] = (1, [gi])

  ll   = [ (k[0], k[1][0]) for k in sorted(cont.items(), key=lambda x: x[1][0], reverse=True) ];
  maxl = max([ len(o[0]) for o in ll ]) + 3 if len(ll) > 0 else 0;

  print "[+] Possible contaminants found in sample %s:" % sn;
  print "  %-*s %-5s" % (maxl, "Organism", "Hits");
  for o in ll:
    print "  %-*s %-5d" % (maxl, o[0], o[1]);
  #efor
  print "\n\n";
  sys.stdout.flush()
  Save(cont, '%s/%s.unmapped_orgs.dat' % (C.outdir, sn))
#efor

sys.exit(0);
