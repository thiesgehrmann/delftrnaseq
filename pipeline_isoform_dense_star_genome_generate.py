#!/usr/bin/python
from pipeline_common import *;

C = init_conf()

cmds = [];

split_gff, split_fasta = cor(C.__isoform_dense_genome_split_output__);

cmds.append("mkdir -p %s" % cor(C.__isoform_dense_genome_generate_dir__));

if C.isoform_dense_build_splice_db:
  cmds.append("STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbFileChrStartEnd %s %s" % (cor(C.C.__isoform_dense_genome_generate_dir__), C.genome, C.__isoform_dense_star_preal_output__(), cor(C.star_gg_opts)));
else:
  cmds.append("STAR --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s %s" % (cor(C.__isoform_dense_genome_generate_dir__), split_fasta, cor(C.star_gg_opts)));
#fi

sys.exit(run_seq_cmds(cmds));

