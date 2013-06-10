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

outputs = zip([ "genes.fpkm_tracking", "isoforms.fpkm_tracking", "transcripts.gtf", "skipped.gtf" ], C.__cufflinks_output__() );
cmds = [];

cmds.append(("cufflinks -o %s " % C.outdir) + \
            ("%s " % C.cufflinks_opts) + \
            (("-G %s " % C.__genome_annot_format_output__()) if C.genome_annot else ("-g %s " % C.genome_guide))  + \
            ("-v ") + \
            ("-b %s " % C.cufflinks_bias_corr if C.cufflinks_bias_corr else "") + \
            ("%s" % C.__pre_cufflinks_sort_output__()) ) ;

for (o , n) in outputs:
  cmds.append("mv '%s/%s' '%s'" % (C.outdir, o, n));
#efor

retval = run_seq_cmds(cmds);
if retval != 0:
  sys.exit(retval);
#fi

sys.exit(retval);

