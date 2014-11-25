#!/usr/bin/env python
from pipeline_common import *;

C = init_conf()

cmds = [];

split_gff = C.outdir + '/split_genome.gff';
out_gff, out_fasta, out_split_info = cor(C.__isoform_dense_genome_split_output__);
star_al_sam        = cor(C.__star_al_output_sam__);
star_al_logs       = cor(C.__star_al_output_logs__);

genome_annot = C.__genome_annot_format_output__();


if C.PE:
  if (cor(C.__isoform_dense_genome_split_is_mean__) == None or cor(C.__isoform_dense_genome_split_is_stdev__) == None):
    print "Getting IS size";
    cmd = "%s/utilities/get_is_dist.sh %s" % (C.inst_loc, ' '.join([ "%s" % sam for sam in star_al_sam ]));
    out = getCommandOutput(cmd);
    IS_mean, IS_stdev = out.split(' ');
  #fi
  if cor(C.__isoform_dense_genome_split_is_mean__) != None:
    IS_mean = cor(C.__isoform_dense_genome_split_is_mean__);
  #fi
  if cor(C.__isoform_dense_genome_split_is_stdev__) != None:
    IS_stdev = str(cor(C.__isoform_dense_genome_split_is_stdev__));
  #fi
else:
  IS_mean  = '0';
  IS_stdev = '0';
#fi



if cor(C.__isoform_dense_genome_split_readlength__) == None:
  print "Getting average readlength"
  cmd = "%s/utilities/get_average_readlength.sh %s" % (C.inst_loc, ' '.join(star_al_logs[2]));
  average_readlength = getCommandOutput(cmd);
else:
  average_readlength = str(cor(C.__isoform_dense_genome_split_readlength__));
#fi


if C.PE:
  cmd = "%s/utilities/split_genome.py split %s %s %s %s %s %s %s %s" % (C.inst_loc, genome_annot, C.genome, average_readlength, IS_mean, IS_stdev, split_gff, out_fasta, out_split_info);
else:
  cmd = "%s/utilities/split_genome.py split %s %s %s 0 0 %s %s" % (C.inst_loc, genome_annot, C.genome, average_readlength, split_gff, out_fasta, out_split_info);
#fi
cmds.append(cmd);

cmds.append("gffread -FE %s -o %s" % (split_gff, out_gff));
cmds.append("rm %s" % split_gff);

print "Splitting genome"
sys.exit(run_seq_cmds(cmds));


