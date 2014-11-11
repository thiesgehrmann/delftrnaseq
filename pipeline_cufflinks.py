#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

outputs = zip([ "genes.fpkm_tracking", "isoforms.fpkm_tracking", "transcripts.gtf", "skipped.gtf" ], C.__cufflinks_output__() );
cmds = [];

cmds.append(("cufflinks -o %s " % C.outdir) + \
            ("%s " % cor(C.cufflinks_opts)) + \
            (("-G %s " % C.__genome_annot_format_output__()) if C.genome_annot else ("-g %s " % C.genome_guide))  + \
            ("-v ") + \
            ("-b %s " % C.cufflinks_bias_corr if C.cufflinks_bias_corr else "") + \
            ("%s" % C.__pre_cufflinks_merge_output__()) ) ;

for (o , n) in outputs:
  cmds.append("mv '%s/%s' '%s'" % (C.outdir, o, n));
#efor

retval = run_seq_cmds(cmds);
if retval != 0:
  sys.exit(retval);
#fi

sys.exit(retval);

