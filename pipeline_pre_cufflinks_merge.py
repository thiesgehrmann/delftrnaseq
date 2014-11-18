#!/usr/bin/env python

from pipeline_common import *;

C = init_conf()

bamin  = C.__star_al_output_bam__();
bamout = C.__pre_cufflinks_merge_output__();

cmd = 'bamtools merge %s -out %s' % (' -in '.join([' '] + bamin), bamout);

sys.exit(run_cmd(cmd));

