#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

prg="STAR"

cmds = [ "mkdir -p %s " % C.__star_pre_splice_output__(),
         "%s --runMode genomeGenerate --genomeDir %s %s --genomeFastaFiles %s" % (prg, C.__star_pre_splice_output__(), cor(C.star_gg_opts), C.genome) ]

sys.exit(run_seq_cmds(cmds));


