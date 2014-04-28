#!/usr/bin/env python
from pipeline_common import *
import os

C = init_conf()

genome_gff = C.__genome_annot_format_output__()
sort_outputs = C.__post_star_al_sort_output__()

cmds = []

for i in xrange(len(sort_outputs)) :
  indexed_bam = sort_outputs[i]
  prefix, ext = os.path.splitext(indexed_bam)
  if C.strand_specific :
    cmd = "htseq-count -i ID -f bam -r pos -s yes -t %s '%s' '%s' >%s.count" % (C.cds_gff_type, indexed_bam, genome_gff, prefix)
  else :
    cmd = "htseq-count -i ID -f bam -r pos -s no -t %s '%s' '%s' >%s.count" % (C.cds_gff_type, indexed_bam, genome_gff, prefix)
  cmds.append(cmd)

print "Covering genome GFF with reads"
sys.stdout.flush()
retval = run_par_cmds(cmds, max_threads = min(C.max_bam_threads, 3), shell = True)

if retval != 0 :
  sys.exit(retval)

sys.exit(0)
