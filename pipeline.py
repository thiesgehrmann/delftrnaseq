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
    fd.write('.PHONY : all\n');
    fd.write('all: %s\n\n' % end)
    fd.write('.PHONY : clean\n');
    fd.write('clean : %s\n' % workdir)
    fd.write('\t@rm -rf %s\n\n' % workdir);


    for (step, script) in self.steps[::-1]:
      #fd.write('.PHONY : %s\n' % step);
      fd.write('%s ${%s_OUT} : ${%s_IN}\n' % (step, step, step));
      fd.write('\t@echo [%s] Starting\n' % step);
      fd.write('\t@${inst_loc}/%s "%s" | tee "${outdir}/%s.std.log"\n' % (script, loc, step));
      fd.write('\t@${inst_loc}/pipeline_notify.py "%s" "%s" "${outdir}/%s.std.log"\n' % (loc, step, step));
      fd.write('\t@echo [%s] Finished\n\n' % step);
    #efor

    fd.close();
  #edef

#eclass

###############################################################################

execfile(os.sys.argv[1]);

C = config()

C.wd_correct()
C.data_check()
C.location = C.pickle();

M = makefile();

M.add_var('inst_loc', C.inst_loc);
M.add_var('outdir', C.outdir);

M.add_step("TRIMMOMATIC", (' '.join(flatten(C.samples))), (' '.join(flatten(C.trimmomatic_output()))), 'pipeline_trimmomatic.py');
M.add_step("STAR_GG", (' '.join(flatten(C.genome))),  C.star_gg_output(), 'pipeline_star_genome_generate.py');

M.add_var("STAR_AL_OUTPUT_NAME", ' '.join(C.star_al_output_name()));
M.add_var("STAR_AL_OUTPUT_CHR", ' '.join(C.star_al_output_chr()));
M.add_var("STAR_AL_OUTPUT_UNMAPPED", ' '.join(flatten(C.star_al_output_unmapped())));
M.add_step("STAR_AL", "${STAR_GG_OUT} ${TRIMMOMATIC_OUT}", "${STAR_AL_OUTPUT_NAME} ${STAR_AL_OUTPUT_CHR} ${STAR_AL_OUTPUT_UNMAPPED}", 'pipeline_star_align.py');

M.add_step("PRE_CUFFLINKS", "${STAR_AL_OUTPUT_NAME}", C.pre_cufflinks_output(), 'pipeline_pre_cufflinks.py');
M.add_step("CUFFLINKS", "${PRE_CUFFLINKS_OUT} %s" % C.genome_annot, ' '.join(C.cufflinks_output()), 'pipeline_cufflinks.py');

M.add_step("CUFFDIFF", "${STAR_AL_OUTPUT_NAME} ${CUFFLINKS_OUT}", ' '.join(C.cuffdiff_output()), 'pipeline_cuffdiff.py');

M.add_step("TRINITY", "${STAR_AL_OUTPUT_UNMAPPED}", ' '.join(C.trinity_output()), 'pipeline_trinity.py');
M.add_step("TRINITY_ORF", "${TRINITY_OUT}", ' '.join(C.trinity_orf_output()), 'pipeline_trinity_orf.py');
M.add_step("UNMAPPED", "${TRINITY_ORF_OUT}", ' '.join(C.unmapped_output()), 'pipeline_unmapped.py');

M.write(C.makefile, C.location, "${CUFFLINKS_OUT}", C.outdir);



