#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

bamin  = C.__post_star_al_output__();
bamout = C.__pre_cufflinks_indiv_sort_output__();
cmds = [];

for i in xrange(len(bamin)):
  bam = bamin[i];
  sn  = C.sample_names[i];
  out = bamout[i][0:-4];

  cmds.append("samtools sort -m 100000000000 '%s' '%s'" % (bam, out));
#efor

retval = run_par_cmds(cmds, max_threads=C.__max_threads__);

if retval != 0:
  sys.exit(retval);
#fi

sys.exit(0);

