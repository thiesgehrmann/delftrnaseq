#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

prg="STAR"

old_gff, old_fasta = cor(C.__isoform_dense_genome_split_output__);
new_gff            = cor(C.__isoform_dense_genome_unsplit_output__);
attr_name          = cor(C.__isoform_dense_attr_name__);
outdir             = cor(C.__isoform_dense_genome_analysis_outdir__);

cmds.append("mkdir -p %s" % outdir);
cmds.append("%s/utlities/splicing_statistics.py %s %s %s %s" % (C.inst_loc, old_gff, new_gff, attr_name, outdir));

sys.exit(run_seq_cmds(cmds));

