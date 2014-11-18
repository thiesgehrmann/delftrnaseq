#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

prg="STAR"

split_gff, split_fasta = cor(C.__isoform_dense_genome_split_output__);

cmds = [ "mkdir -p %s " % cor(C.__isoform_dense_star_pre_splice_output__),
         "%s --runMode genomeGenerate --genomeDir %s %s --genomeFastaFiles %s" % (prg, cor(C.__isoform_dense_star_pre_splice_output__), cor(C.star_gg_opts), split_fasta) ]

sys.exit(run_seq_cmds(cmds));


