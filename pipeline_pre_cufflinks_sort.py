#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

bamin  = C.__pre_cufflinks_merge_output__();
bamout = C.__pre_cufflinks_sort_output__()[0:-4];

print bamin, bamout;
sys.stdout.flush();

cmd = "samtools sort -m 100000000000 '%s' '%s'" % (bamin, bamout);

sys.exit(run_cmd(cmd));

