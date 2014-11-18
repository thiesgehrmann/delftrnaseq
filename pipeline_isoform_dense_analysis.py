#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

prg="STAR"

cmds = [ "mkdir -p %s " % C.__star_gg_output__() ];

if C.build_splice_db:
  cmds.append("%s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbFileChrStartEnd %s %s" % (prg, C.__star_gg_output__(), C.genome, C.__star_preal_output__(), cor(C.star_gg_opts)));
else:
  cmds.append("%s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s %s" % (prg, C.__star_gg_output__(), C.genome, cor(C.star_gg_opts)));
#fi

sys.exit(run_seq_cmds(cmds));


