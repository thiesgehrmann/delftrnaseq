#!/usr/bin/env python
from pipeline_common import *;
C =  init_conf()

outputs = [ "genes.fpkm_tracking", "isoforms.fpkm_tracking", "transcripts.gtf", "skipped.gtf" ];
rnames  = C.__cufflinks_indiv_output__();
bamin = C.__post_star_al_sort_output__();

cmds = [];

for i in xrange(len(bamin)):

  cmds.append(("cufflinks -o %s " % C.outdir) + \
              ("%s " % cor(C.cufflinks_indiv_opts)) + \
              (("-G %s " % C.__genome_annot_format_output__()) if C.genome_annot else ("-g %s " % C.genome_guide))  + \
              ("-v ") + \
              ("-b %s " % C.cufflinks_bias_corr if C.cufflinks_bias_corr else "") + \
              ("%s" % bamin[i]) ) ;


  outfs = zip(outputs, rnames[i]);

  for (o , n) in outfs:
    cmds.append("mv '%s/%s' '%s'" % (C.outdir, o, n));
  #efor

#efor


retval = run_seq_cmds(cmds);

if retval != 0:
  sys.exit(retval);
#fi

sys.exit(retval);
