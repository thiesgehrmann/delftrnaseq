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

bf  = C.__post_star_al_sort_output__();
gtf = C.__cufflinks_output__()[2];
out = C.__cuffdiff_output__();
cuff_files = [ "bias_params.info", "cds.count_tracking", "cds.diff", "cds.fpkm_tracking", "cds.read_group_tracking", "cds_exp.diff", "gene_exp.diff", "genes.count_tracking", "genes.fpkm_tracking", "genes.read_group_tracking", "isoform_exp.diff", "isoforms.count_tracking", "isoforms.fpkm_tracking", "isoforms.read_group_tracking", "promoters.diff", "read_groups.info", "run.info", "splicing.diff", "tss_group_exp.diff", "tss_groups.count_tracking", "tss_groups.fpkm_tracking", "tss_groups.read_group_tracking", "var_model.info" ];


cmds = [];

for i in xrange(len(C.cuffdiff_cmp):
  a, b = C.cuffdiff_cmp[i];
  a_reps = [ bf[j] for j in xrange(len(bf)) if C.sample_labels[j] == a ];
  b_reps = [ bf[j] for j in xrange(len(bf)) if C.sample_labels[j] == b ];
  files = out[i];
  

  cmd = "cuffdiff " \
      + "%s " % cor(C.cuffdiff_opts) \
      + "-L '%s,%s' " % (str(a), str(b)) \
      + "-b '%s' " % C.genome \
      + "%s " % gtf \
      + "%s " % ','.join(a_reps) \
      + "%s " % ','.join(b_reps);

  cmds.append(cmd);

  for i in xrange(len(cuff_files)):
    cmd = "mv %s %s" % (cuff_files[i], files[i]);
    cmds.append(cmd); 


#efor

#retval = run_seq_cmds(cmds);

print cmds[0];


sys.exit(retval);
