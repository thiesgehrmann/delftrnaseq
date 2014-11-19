#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

genome_annot = C.__genome_annot_format_output__();

split_gff, split_fasta = cor(C.__isoform_dense_genome_split_output__);

cuff_outputs = cor(C.__isoform_dense_cufflinks_output__);

unsplit_gff = "%s/unsplit_genome.gff" % C.outdir;

out_gff     = cor(C.__isoform_dense_genome_unsplit_output__);

cmds = [];

cmds.append("%s/utilities/split_genome.py unsplit %s %s %s %s %s" % (C.inst_loc, genome_annot, split_fasta, cuff_outputs[-1], unsplit_gff));
cmds.append("gffread -EF %s -o %s" % (unsplit_gff, out_gff));
cmds.append("rm %s" % unsplit_gff);

sys.exit(run_seq_cmds(cmds));


