#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

SORT = C.__star_al_output_bam__();

cmds = []

for i in xrange(len(SORT)):
  s  = SORT[i][:-4];
  cmds.append("samtools index '%s.bam'" % s);

print "Indexing BAM files"; sys.stdout.flush();
retval = run_par_cmds(cmds, max_threads = min(C.max_bam_threads, 2));
if retval != 0 :
  sys.exit(retval);

sys.exit(0);

