#!/usr/bin/python

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
    for (step, script) in self.steps:
      fd.write('\t@printf "%%-20s" "[%s]: "; ' % step);
      fd.write('ls ${%s_OUT} &>/dev/null; ' % step);
      fd.write('if [ $$? -eq 0 ]; then');
      fd.write('  echo "COMPLETE";');
      fd.write('elif [ -n "`ps aux | grep \'${inst_loc}/%s %s\' | grep -v \'grep\'`" ]; then' % (script, loc));
      fd.write('  c=`ls ${%s_OUT} 2> /dev/null | wc -l`;' % step);
      fd.write('  t=`echo ${%s_OUT} | tr \' \' \'\\n\' | wc -l`;' % step);
      fd.write('  echo "RUNNING $$c / $$t";');
      fd.write('else');
      fd.write('  echo "INCOMPLETE";');
      fd.write('fi\n');
    #efor
    fd.write('.PHONY : touch\n');
    fd.write('touch: status\n');
    for (step, script) in self.steps:
      fd.write('\t@ls ${%s_OUT} &>/dev/null; if [ $$? -eq 0 ]; then touch ${%s_OUT}; fi\n' % (step, step));
    #efor
    fd.write('\n\n');


    for (step, script) in self.steps[::-1]:
      #fd.write('.PHONY : %s\n' % step);
      fd.write('%s ${%s_OUT} : ${%s_IN}\n' % (step, step, step));
      fd.write('\t@echo [%s] Starting\n' % step);
      fd.write('\t@${inst_loc}/%s "%s" | tee "${outdir}/%s.std.log";' % (script, loc, step));
      fd.write(' ps=$${PIPESTATUS[0]}; if [ $$ps -eq 0 ]; then stat="COMPLETED"; else stat="FAILED"; fi;');
      fd.write(' ${inst_loc}/pipeline_notify.py "%s" "%s" "${outdir}/%s.std.log" "$$stat"; echo "[%s] $$stat" ; exit "$$ps";\n' % (loc, step, step, step));
    #efor

    fd.close();
  #edef

#eclass

###############################################################################

execfile(os.sys.argv[1]);

C = config()
C.inst_loc = install_loc;
C.wd_correct()
C.data_check()
C.location = C.pickle();

M = makefile();

M.add_var('inst_loc', C.inst_loc);
M.add_var('outdir', C.outdir);

M.add_step("TRIMMOMATIC", (' '.join(flatten(C.samples))), (' '.join(flatten(C.__trimmomatic_output__()))), 'pipeline_trimmomatic.py');
M.add_step("STAR_GG", (' '.join(flatten(C.genome))),  C.__star_gg_output__(), 'pipeline_star_genome_generate.py');

M.add_var("STAR_AL_OUTPUT_SAM", ' '.join(C.__star_al_output_sam__()));
M.add_var("STAR_AL_OUTPUT_UNMAPPED", ' '.join(flatten(C.__star_al_output_unmapped__())));
M.add_step("STAR_AL", "${STAR_GG_OUT} ${TRIMMOMATIC_OUT}", "${STAR_AL_OUTPUT_SAM} ${STAR_AL_OUTPUT_UNMAPPED}", 'pipeline_star_align.py');
M.add_step("POST_STAR_AL", "${STAR_AL_OUTPUT_SAM}", ' '.join(C.__post_star_al_output__()), 'pipeline_post_star_al.py');

M.add_step("PRE_CUFFLINKS", "${POST_STAR_AL_OUT}", C.__pre_cufflinks_output__(), 'pipeline_pre_cufflinks.py');
M.add_step("CUFFLINKS", "${PRE_CUFFLINKS_OUT} %s" % C.genome_annot, ' '.join(C.__cufflinks_output__()), 'pipeline_cufflinks.py');

M.add_step("CUFFDIFF", "${POST_STAR_AL_OUT} ${CUFFLINKS_OUT}", ' '.join(C.__cuffdiff_output__()), 'pipeline_cuffdiff.py');

M.add_step("TRINITY", "${STAR_AL_OUTPUT_UNMAPPED}", ' '.join(C.__trinity_output__()), 'pipeline_trinity.py');
M.add_step("TRINITY_ORF", "${TRINITY_OUT}", ' '.join(C.__trinity_orf_output__()), 'pipeline_trinity_orf.py');
M.add_step("UNMAPPED_BLAST", "${TRINITY_ORF_OUT} %s" % (C.blast_db if C.blast_db != None else ""), ' '.join(C.__unmapped_blast_output__()), 'pipeline_unmapped_blast.py');
M.add_step("UNMAPPED", "${UNMAPPED_BLAST}", ' '.join(C.__unmapped_output__()), 'pipeline_unmapped.py');

M.write(C.makefile, C.location, "${CUFFLINKS_OUT}", C.outdir);



