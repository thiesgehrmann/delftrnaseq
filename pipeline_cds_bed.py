#!/usr/bin/env python
from pipeline_common import *
import os

C = init_conf()

genome_gff = C.__genome_annot_format_output__()
genome_bed = C.__cds_bed_genome_output__()
sort_outputs = C.__post_star_al_sort_output__()

cmd = "gff2bed <%s | grep 'CDS'>%s" % (genome_gff, genome_bed)
print "Creating genome BED"
sys.stdout.flush()

retval = run_shell(cmd)
if retval != 0 :
  sys.exit(retval)

cmds = []

for i in xrange(len(sort_outputs)) :
  indexed_bam = sort_outputs[i]
  prefix, ext = os.path.splitext(indexed_bam)
  cmd = "multiBamCov -p -bams '%s' -bed '%s' >%s.bed" % (indexed_bam, genome_bed, prefix)
  cmds.append(cmd)
  #retval = run_shell(cmd)
  #if retval != 0 :
  #    sys.exit(retval)

print "Covering genome BED with reads"
sys.stdout.flush()
retval = run_par_cmds(cmds, max_threads = min(C.max_bam_threads, 2), shell = True)

if retval != 0 :
  sys.exit(retval)

sys.exit(0)
