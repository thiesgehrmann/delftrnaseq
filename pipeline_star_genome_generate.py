#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

prg="STAR"

cmds = [ "mkdir -p %s " % C.__star_gg_output__(),
         "%s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbFileChrStartEnd %s %s" % (prg, C.__star_gg_output__(), C.genome, C.__star_preal_output__(), cor(C.star_gg_opts))]

sys.exit(run_seq_cmds(cmds));


