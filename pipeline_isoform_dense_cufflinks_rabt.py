#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

cmds = [];

BAM = cor(C.__isoform_dense_star_align_output_merged__);
split_gff, split_fasta, split_info = cor(C.__isoform_dense_genome_split_output__);

cuff_outdir  = cor(C.__isoform_dense_cufflinks_outdir__)
cuff_outputs = cor(C.__isoform_dense_cufflinks_output__);
cuffopts = cor(C.isoform_dense_cufflinks_opts);

attr_name = cor(C.__isoform_dense_attr_name__);

cmds.append("mkdir -p %s" % cuff_outdir);

cmds.append(("cufflinks -o %s " % cuff_outdir) + \
            ("%s " % cuffopts ) + \
            ("-g %s " % split_gff) +\
            ("-v ") + \
            ("-b %s " % C.cufflinks_bias_corr if C.cufflinks_bias_corr else "") + \
            ("%s" % BAM) ) ;

cmds.append("gffread -EF %s -o %s" % (cuff_outputs[2], cuff_outputs[-2]));

tmp_gff = cuff_outputs[-1]+ '.tmp.gff';

cmds.append("%s/utilities/split_genome.py setx2y %s %s attribute seqname name %s" % (C.inst_loc, cuff_outputs[-2], "geneID", tmp_gff));

cmds.append("%s/utilities/split_genome.py setx2y %s %s attribute seqname name %s" % (C.inst_loc, tmp_gff, attr_name, cuff_outputs[-1]));

cmds.append("rm %s" % tmp_gff);

sys.exit(run_seq_cmds(cmds));


