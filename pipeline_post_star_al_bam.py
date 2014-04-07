#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

ALN = C.__star_al_output_sam__();
BAM = C.__post_star_al_bam_output__();

cmds = [];

for i in xrange(len(ALN)):
  a = ALN[i];
  b = BAM[i];

  cmds.append("samtools view -Sb '%s' -o '%s'" % (a, b));

#efor

print "Converting SAM to BAM"; sys.stdout.flush();
retval = run_par_cmds(cmds, max_threads=C.max_bam_threads);
if retval != 0:
 sys.exit(retval);
#fi 

sys.exit(0);

