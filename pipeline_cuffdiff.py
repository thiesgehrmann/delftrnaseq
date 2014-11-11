#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

bf  = C.__star_al_output_bam__();
gtf = C.__cufflinks_output__()[2];
out = C.__cuffdiff_output__();
cuff_files = [ "bias_params.info", "cds.count_tracking", "cds.diff", "cds.fpkm_tracking", "cds.read_group_tracking", "cds_exp.diff", "gene_exp.diff", "genes.count_tracking", "genes.fpkm_tracking", "genes.read_group_tracking", "isoform_exp.diff", "isoforms.count_tracking", "isoforms.fpkm_tracking", "isoforms.read_group_tracking", "promoters.diff", "read_groups.info", "run.info", "splicing.diff", "tss_group_exp.diff", "tss_groups.count_tracking", "tss_groups.fpkm_tracking", "tss_groups.read_group_tracking", "var_model.info" ];

cmds = [];

all_reps = []
for i in xrange(len(C.label_names)):
    all_reps.append([ bf[j] for j in xrange(len(bf)) if C.sample_labels[j] == i ]);

labels = ",".join([str(i) for i in range(len(C.label_names))])

cmd = "cuffdiff " \
    + "%s " % cor(C.cuffdiff_opts) \
    + "-L '%s' " % labels \
    + "-b '%s' " % C.genome \
    + "%s " % gtf

for reps in all_reps:
    cmd = cmd + "%s " % ','.join(reps) 

cmds.append(cmd);

for i in xrange(len(cuff_files)):
    cmd = "mv %s %s" % (cuff_files[i], out[i]);
    cmds.append(cmd); 
   

retval = run_seq_cmds(cmds);

sys.exit(retval);
