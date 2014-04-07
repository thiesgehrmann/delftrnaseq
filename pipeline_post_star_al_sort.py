#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

BAM  = C.__post_star_al_bam_output__();
SORT = C.__post_star_al_sort_output__();

cmds = []

for i in xrange(len(BAM)):
  b  = BAM[i];
  s  = SORT[i][:-4];

  cmds.append("samtools sort -m 50000000000 '%s' '%s'" % (b, s));

#efor

print "Sorting BAM files"; sys.stdout.flush();
retval = run_par_cmds(cmds, max_threads=min(C.max_bam_threads,2));
if retval != 0:
  sys.exit(retval);
#fi

sys.exit(0);

