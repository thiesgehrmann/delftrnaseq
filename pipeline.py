#!/usr/bin/env python

import os
import inspect
import tempfile
import pickle

###############################################################################

install_loc = os.path.realpath(os.path.dirname(os.sys.argv[0]));
os.sys.path.append(install_loc);

from pipeline_common import *

###############################################################################
###############################################################################
###############################################################################

class makefile:

  vars = [];
  steps = []

  def add_step(self, name, IN, OUT, script):
    self.add_var('%s_IN' % name, IN);
    self.add_var('%s_OUT' % name, OUT);
    self.steps = self.steps + [ ( name, script ) ];
  #edef

  def add_var(self, name, val):
    self.vars = self.vars + [ ( name, val ) ];
  #edef

  def write(self, mfile, loc, end, workdir):
    fd = open(mfile, 'w');

    for (name, val) in self.vars:
      fd.write('%s = %s\n' % (name, val));
    #efor

    fd.write('\n\n');
    fd.write('all: %s\n\n' % end)
    fd.write('.PHONY : clean\n');
    fd.write('clean : %s\n' % workdir)
    fd.write('\t@rm -rf %s\n' % workdir);

    fd.write('.PHONY : status\n');
    fd.write('status: \n');
    maxl = max([len(n) + 4 for (n,s) in self.steps]);
    for (step, script) in self.steps:
      fd.write('\t@printf "%%-*s" %d "[%s]: "; ' % (maxl, step));
      fd.write('ls ${%s_OUT} &>/dev/null; ' % step);
      fd.write('res=$$?; ');
      fd.write('if [ -n "`ps aux | grep \'${inst_loc}/%s %s\' | grep -v \'grep\'`" ]; then' % (script, loc));
      fd.write('  c=`ls ${%s_OUT} 2> /dev/null | wc -l`;' % step);
      fd.write('  t=`echo ${%s_OUT} | tr \' \' \'\\n\' | wc -l`;' % step);
      fd.write('  printf "RUNNING %3s / %3s\\n" "$$c" "$$t";');
      fd.write('elif [ $$res -eq 0 ]; then');
      fd.write('  echo "COMPLETE";');
      fd.write('else');
      fd.write('  echo "INCOMPLETE";');
      fd.write('fi\n');
    #efor
    fd.write('.PHONY : touch\n');
    fd.write('touch: status\n');
    fd.write('\t@find ${outdir} | xargs touch\n')
    fd.write('\t@find . | xargs touch\n')
    for (step, script) in self.steps:
      fd.write('\t@ls ${%s_OUT} &>/dev/null;' % step);
      fd.write('if [ $$? -eq 0 ]; then ');
      fd.write('  touch ${%s_OUT};' % step);
      fd.write('fi\n');
    #efor
    fd.write('\n\n');

    for (step, script) in self.steps[::-1]:
      fd.write('%s ${%s_OUT} : ${%s_IN}\n' % (step, step, step));
      fd.write('\t@echo [%s] Starting at `date`\n' % step);
      fd.write('\t@rm -rf ${%s_OUT}\n' % step);
      fd.write('\t@ds=`date`; ');
      fd.write('${inst_loc}/%s "%s" 2>&1 | tee "${outdir}/%s.std.log";' % (script, loc, step));
      fd.write('ps=$${PIPESTATUS[0]}; ');
      fd.write('if [ $$ps -eq 0 ]; then ');
      fd.write('  stat="COMPLETED"; ');
      fd.write('else ');
      fd.write('  stat="FAILED"; ');
      fd.write('fi; ');
      fd.write('mkdir -p ${outdir}; ');
      fd.write('${inst_loc}/pipeline_notify.py "%s" "%s" "${outdir}/%s.std.log" "$$stat" "$$ds" "`date`"; echo "[%s] $$stat at `date`" ; exit $$ps;\n' % (loc, step, step, step));
    #efor

    fd.close();
  #edef

#eclass

###############################################################################

execfile(os.sys.argv[1]);

try:
  C = config();
except:
  print "ERROR: You must define a class 'config(PIPELINECONF)' in your '%s' file." % os.sys.argv[1];
  sys.exit(1);
#etry

if C.inst_loc == None:
  C.inst_loc = install_loc;
#fi

C.wd_correct();
C.data_check();
C.location = C.pickle();

###############################################################################

M = makefile();

M.add_var('inst_loc', C.inst_loc);
M.add_var('outdir', C.outdir);

  # Quality step
M.add_step("TRIMMOMATIC", (' '.join(flatten(C.samples))), (' '.join(flatten(C.__trimmomatic_output__()))), 'pipeline_trimmomatic.py');

  # Genome annotation step
if (C.genome_annot == None) and (C.genome_guide == None):
  M.add_step("GENOME_ANNOT", "${TRIMMOMATIC} %s" % C.genome, C.__genome_annot_output__(), 'pipeline_genome_annot.py');
#fi

  # STAR alignment steps
if C.build_splice_db:
  M.add_step("STAR_PRE_SPLICE", C.genome, C.__star_pre_splice_output__(), 'pipeline_star_pre_splice.py');
  M.add_step("STAR_SPLICE", "${STAR_PRE_SPLICE_OUT} ${TRIMMOMATIC_OUT}", C.__star_preal_output__(), 'pipeline_star_splice.py');
  M.add_step("STAR_GG", "${STAR_SPLICE_OUT}",  C.__star_gg_output__(), 'pipeline_star_genome_generate.py');
else:
  M.add_step("STAR_GG", C.genome, C.__star_gg_output__(), 'pipeline_star_genome_generate.py');
#fi

M.add_var("STAR_AL_OUTPUT_SAM", ' '.join(C.__star_al_output_sam__()));
M.add_var("STAR_AL_OUTPUT_UNMAPPED", ' '.join(flatten(C.__star_al_output_unmapped__())));
M.add_step("STAR_AL", "${STAR_GG_OUT} ${TRIMMOMATIC_OUT}", "${STAR_AL_OUTPUT_SAM} ${STAR_AL_OUTPUT_UNMAPPED}", 'pipeline_star_align.py');
M.add_step("POST_STAR_AL_BAM", "${STAR_AL_OUTPUT_SAM}", ' '.join(C.__post_star_al_bam_output__()), 'pipeline_post_star_al_bam.py');
M.add_step("POST_STAR_AL_SORT", "${POST_STAR_AL_BAM_OUT}", ' '.join(C.__post_star_al_sort_output__()), 'pipeline_post_star_al_sort.py');

  # CUFF steps
M.add_step("PRE_CUFFLINKS_MERGE", "${POST_STAR_AL_BAM_OUT}", C.__pre_cufflinks_merge_output__(), 'pipeline_pre_cufflinks_merge.py');
M.add_step("PRE_CUFFLINKS_SORT", "${PRE_CUFFLINKS_MERGE_OUT}", C.__pre_cufflinks_sort_output__(), 'pipeline_pre_cufflinks_sort.py');
M.add_step("GENOME_ANNOT_FORMAT", C.genome_annot, C.__genome_annot_format_output__(), 'pipeline_genome_annot_format.py');
M.add_step("CUFFLINKS", "${PRE_CUFFLINKS_SORT_OUT} ${GENOME_ANNOT_FORMAT_OUT}", ' '.join(C.__cufflinks_output__()), 'pipeline_cufflinks.py');
M.add_step("CUFFDIFF", "${POST_STAR_AL_SORT_OUT} ${CUFFLINKS_OUT}", ' '.join(flatten(C.__cuffdiff_output__())), 'pipeline_cuffdiff.py');

  # CUFF_INDIV steps
#M.add_step("PRE_CUFFLINKS_INDIV_SORT", "${POST_STAR_AL_BAM_OUT}", ' '.join(C.__pre_cufflinks_indiv_sort_output__()), 'pipeline_pre_cufflinks_indiv_sort.py');
M.add_step("CUFFLINKS_INDIV", "${POST_STAR_AL_SORT_OUT} ${GENOME_ANNOT_FORMAT_OUT}", ' '.join(flatten(C.__cufflinks_indiv_output__())), 'pipeline_cufflinks_indiv.py');

  # Contamination steps
M.add_step("TRINITY", "${STAR_AL_OUTPUT_UNMAPPED}", ' '.join(C.__trinity_output__()), 'pipeline_trinity.py');
M.add_step("TRINITY_ORF", "${TRINITY_OUT}", ' '.join(C.__trinity_orf_output__()), 'pipeline_trinity_orf.py');
M.add_step("UNMAPPED_BLAST", "${TRINITY_ORF_OUT} %s" % (C.blast_db if C.blast_db != None else ""), ' '.join(C.__unmapped_blast_output__()), 'pipeline_unmapped_blast.py');
M.add_step("UNMAPPED", "${UNMAPPED_BLAST_OUT}", ' '.join(C.__unmapped_output__()), 'pipeline_unmapped.py');
M.add_step("CUFFDIFF_COMBINE", "${CUFFDIFF_OUT}", ' '.join(C.__cuffdiff_combine_output__()), 'pipeline_cuffdiff_combine.py');

  # REPORTS
M.add_step("QUALITYREPORT", "${TRIMMOMATIC_OUT} ${STAR_AL_OUT}", ' '.join(flatten(C.__quality_output__())), 'pipeline_quality.py');
M.add_step("POSTANALYSIS", "${CUFFDIFF_COMBINE_OUT}", ' '.join(flatten(C.__analysis_output__())), 'pipeline_analysis.py');

M.write(C.makefile, C.location, "${POSTANALYSIS}", C.outdir);

###############################################################################

sys.exit(0);
