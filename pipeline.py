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
    self.add_var('%s_LOG' % name, "${outdir}/%s.std.log" % name);
    self.add_var('%s_OUT' % name, OUT + ' ${%s_LOG}' % name);
    self.steps = self.steps + [ ( name, script ) ];
  #edef

  def add_var(self, name, val):
    self.vars = self.vars + [ ( name, val ) ];
  #edef

  def write(self, mfile, loc, end, workdir, title):
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
    status_header = "PIPELINE STATUS FOR %s" % title;
    spacing = ' ' * (((maxl + 19) - len(status_header)) / 2);
    fd.write('\t@echo "%s%s";' % (spacing, status_header));
    fd.write('echo "%s%s"; echo "";\n' % (spacing, '-' * len(status_header)));
    for (step, script) in self.steps:
      fd.write('\t@printf "%%-*s" %d "[%s]: "; ' % (maxl, step));
      fd.write('ls ${%s_OUT} &>/dev/null; ' % step);
      fd.write('complete_res=$$?; ');
      fd.write('ls ${%s_IN} &>/dev/null;' % step);
      fd.write('ready_res=$$?;');
      fd.write('in_newest=`for x in ${%s_IN}; do if [ -e "$$x" ]; then echo "$$x"; fi; done | xargs stat -c \'%%Y\' 2> /dev/null | sort | tail -n1`;' % step);
      fd.write('for x in ${%s_IN}; do if [ -e "$$x" ]; then echo "$$x"; fi; done | xargs stat -c \'%%Y\' 2> /dev/null | sort | tail -n1 &> /dev/null;' % step);
      fd.write('in_newest_status=$${PIPESTATUS[1]};');
      fd.write('if [ ! $$in_newest_status -eq 0 ]; then');
      fd.write('  in_newest=`date \'+%s\'`;');
      fd.write('fi;');
      fd.write('out_oldest=`for x in ${%s_OUT}; do if [ -e "$$x" ]; then echo $$x; fi; done | xargs stat -c \'%%Y\' 2> /dev/null | sort | head -n1`;' % step);
      fd.write('for x in ${%s_OUT}; do if [ -e "$$x" ]; then echo $$x; fi; done | xargs stat -c \'%%Y\' 2> /dev/null | sort | tail -n1 &> /dev/null;' % step);
      fd.write('out_oldest_status=$${PIPESTATUS[1]};');
      fd.write('if [ ! $$out_oldest_status -eq 0 ]; then');
      fd.write('  out_oldest=`date \'+%s\'`;');
      fd.write('fi;');
      fd.write('if [ "$$ready_res" -eq 0 ]; then');
      fd.write('  echo -en "READY   ";');
      fd.write('else');
      fd.write('  echo -en "UNREADY ";');
      fd.write('fi;');
      fd.write('if [ -n "`ps aux | grep \'${inst_loc}/%s %s\' | grep -v \'grep\'`" ]; then' % (script, loc));
      fd.write('  c=`ls -d ${%s_OUT} 2> /dev/null | wc -l`;' % step);
      fd.write('  t=`echo ${%s_OUT} | tr \' \' \'\\n\' | wc -l`;' % step);
      fd.write('  printf "RUNNING %3s / %3s\\n" "$$c" "$$t";');
      fd.write('elif [ $$complete_res -eq 0 ]; then');
      fd.write('  if [ $$in_newest -gt $$out_oldest ]; then');
      fd.write('    echo "REPEAT";');
      fd.write('  else');
      fd.write('    echo "COMPLETE";');
      fd.write('  fi;');
      fd.write('else');
      fd.write('  echo "INCOMPLETE";');
      fd.write('fi;');
      #fd.write('echo -$$out_oldest-, -$$in_newest-;');
      fd.write('\n');
    #efor
    fd.write('\t@echo "";');
    fd.write('echo "READY      - The direct dependencies of this target have been met.";');
    fd.write('echo "UNREADY    - The direct dependencies of this target have not been met.";');
    fd.write('echo "COMPLETE   - The products of this target have been created more recently than its direct dependencies";');
    fd.write('echo "REPEAT     - The products of this target have been created, but before its direct dependencies. You may want to run touch to resolve this.";');
    fd.write('echo "RUNNING    - This target is currently running.";');
    fd.write('echo "INCOMPLETE - The products of this target have not (all) been created."\n');

    fd.write('.PHONY : just_touch\n');
    fd.write('just_touch:\n');
    fd.write('\t@find ${outdir} -print0 | xargs -0r touch;\n')
    #fd.write('\tfind . | xargs touch -a\n')
    for (step, script) in self.steps:
      fd.write('\t@ls ${%s_IN} &>/dev/null;' % step);
      fd.write('if [ $$? -eq 0 ]; then ');
      fd.write('  echo -e ${%s_IN} | xargs touch;' % step);
      fd.write('fi\n');
    #efor
    for (step, script) in self.steps:
      fd.write('\t@ls ${%s_OUT} &>/dev/null;' % step);
      fd.write('if [ $$? -eq 0 ]; then ');
      fd.write('  echo -e ${%s_OUT} | xargs touch;' % step);
      fd.write('fi\n');
    #efor

    fd.write('.PHONY : touch\n');
    fd.write('touch: just_touch status\n');

    for (step, script) in self.steps:
      fd.write('clean_%s: status\n' % step);
      fd.write('\t@rm -rf ${%s_OUT} &>/dev/null;\n' % step);
    #efor
    fd.write('\n\n');

    for (step, script) in self.steps[::-1]:
      fd.write('%s ${%s_OUT} : ${%s_IN}\n' % (step, step, step));
      fd.write('\t@echo [%s] Starting at `date`\n' % step);
      fd.write('\t@rm -rf ${%s_OUT}\n' % step);
      fd.write('\t@ds=`date`; ');
      fd.write('${inst_loc}/%s "%s" 2>&1 | tee "${%s_LOG}";' % (script, loc, step));
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

M.add_var("STAR_AL_OUTPUT_BAM", ' '.join(C.__star_al_output_bam__()));
M.add_var("STAR_AL_OUTPUT_UNMAPPED", ' '.join(flatten(C.__star_al_output_unmapped__())));
if C.isoform_dense_analysis:
  M.add_var("STAR_AL_OUTPUT_SAM", ' '.join(C.__star_al_output_sam__()));
  M.add_step("STAR_AL", "${STAR_GG_OUT} ${TRIMMOMATIC_OUT}", "${STAR_AL_OUTPUT_BAM} ${STAR_AL_OUTPUT_SAM} ${STAR_AL_OUTPUT_UNMAPPED}", 'pipeline_star_align.py');
else:
  M.add_step("STAR_AL", "${STAR_GG_OUT} ${TRIMMOMATIC_OUT}", "${STAR_AL_OUTPUT_BAM} ${STAR_AL_OUTPUT_UNMAPPED}", 'pipeline_star_align.py');
#fi
#M.add_step("POST_STAR_AL_BAM", "${STAR_AL_OUTPUT_SAM}", ' '.join(C.__post_star_al_bam_output__()), 'pipeline_post_star_al_bam.py');
#M.add_step("POST_STAR_AL_SORT", "${POST_STAR_AL_BAM_OUT}", ' '.join(C.__post_star_al_sort_output__()), 'pipeline_post_star_al_sort.py');

  # CUFF steps
M.add_step("PRE_CUFFLINKS_MERGE", "${STAR_AL_OUTPUT_BAM}", C.__pre_cufflinks_merge_output__(), 'pipeline_pre_cufflinks_merge.py');
#M.add_step("PRE_CUFFLINKS_SORT", "${PRE_CUFFLINKS_MERGE_OUT}", C.__pre_cufflinks_sort_output__(), 'pipeline_pre_cufflinks_sort.py');
M.add_step("GENOME_ANNOT_FORMAT", C.genome_annot, C.__genome_annot_format_output__(), 'pipeline_genome_annot_format.py');
M.add_step("CUFFLINKS", "${PRE_CUFFLINKS_MERGE_OUT} ${GENOME_ANNOT_FORMAT_OUT}", ' '.join(C.__cufflinks_output__()), 'pipeline_cufflinks.py');
M.add_step("CUFFDIFF", "${STAR_AL_OUTPUT_BAM} ${CUFFLINKS_OUT}", ' '.join(flatten(C.__cuffdiff_output__())), 'pipeline_cuffdiff.py');

  # Contamination steps
if C.check_contamination :
  M.add_var("UNMAPPED_OUT", ' '.join(C.__unmapped_output__()));
  M.add_step("TRINITY", "${STAR_AL_OUTPUT_UNMAPPED}", ' '.join(C.__trinity_output__()), 'pipeline_trinity.py');
  M.add_step("TRINITY_ORF", "${TRINITY_OUT}", ' '.join(C.__trinity_orf_output__()), 'pipeline_trinity_orf.py');
  M.add_step("UNMAPPED_BLAST", "${TRINITY_ORF_OUT} %s" % (C.blast_db if C.blast_db != None else ""), ' '.join(C.__unmapped_blast_output__()), 'pipeline_unmapped_blast.py');
  M.add_step("UNMAPPED", "${UNMAPPED_BLAST_OUT}", ' '.join(C.__unmapped_output__()), 'pipeline_unmapped.py');
  #M.add_step("CONTAMINATIONREPORT", "${UNMAPPED_OUT}" '', '');
#fi

  # REPORTS
if C.perform_quality_report:

  M.add_step("POST_STAR_AL_INDEX", "${STAR_AL_OUTPUT_BAM}", ' '.join(C.__post_star_al_index_output__()), 'pipeline_post_star_al_index.py');
  M.add_step("CDS_GFF", "${POST_STAR_AL_INDEX_OUT} ${GENOME_ANNOT_FORMAT_OUT}", ' '.join(C.__cds_gff_output__()), 'pipeline_cds_gff.py');
  M.add_step("READ_DISTRIBUTION", "${STAR_AL_OUTPUT_BAM}", ' '.join(C.__read_distribution_output__()), 'pipeline_read_distribution.py');
  M.add_step("FASTQC", "${STAR_AL_OUT} ${STAR_AL_OUTPUT_BAM}", ' '.join(C.__fastqc_output__()), 'pipeline_fastqc.py');

  if C.check_contamination :
    M.add_step("QUALITYREPORT", "${TRIMMOMATIC_OUT} ${STAR_AL_OUT} ${FASTQC_OUT} ${CDS_GFF_OUT} ${READ_DISTRIBUTION_OUT} ${UNMAPPED_OUT}", ' '.join(flatten(C.__quality_output__())), 'pipeline_quality.py');
  else :
    M.add_step("QUALITYREPORT", "${TRIMMOMATIC_OUT} ${STAR_AL_OUT} ${FASTQC_OUT} ${CDS_GFF_OUT} ${READ_DISTRIBUTION_OUT}", ' '.join(flatten(C.__quality_output__())), 'pipeline_quality.py');
  #fi
#fi

if C.perform_analysis:
  M.add_var('CUFFDIFF_COMBINE_ANNOTATION_FILES', ' '.join(C.annotation_files));
  M.add_step("CUFFDIFF_COMBINE", "${CUFFDIFF_OUT} ${CUFFDIFF_COMBINE_ANNOTATION_FILES}", ' '.join(C.__cuffdiff_combine_output__()), 'pipeline_cuffdiff_combine.py');
  M.add_step("POSTANALYSIS", "${CUFFDIFF_COMBINE_OUT}", ' '.join(flatten(C.__analysis_output__())), 'pipeline_analysis.py');
#fi

if C.isoform_dense_analysis:
  M.add_step("ISOFORM_DENSE_SPLIT_GENOME", "${STAR_AL_OUTPUT_SAM} ${GENOME_ANNOT_FORMAT_OUT}", ' '.join(C.__isoform_dense_genome_split_output__()), 'pipeline_isoform_dense_genome_split.py');
  if C.isoform_dense_build_splice_db:
    M.add_step("ISOFORM_DENSE_STAR_PRE_SPLICE", "${ISOFORM_DENSE_SPLIT_GENOME_OUT}", cor(C.__isoform_dense_star_pre_splice_output__), 'pipeline_isoform_dense_star_pre_splice.py');
    M.add_step("ISOFORM_DENSE_STAR_SPLICE", "${ISOFORM_DENSE_STAR_PRE_SPLICE_OUT} ${TRIMMOMATIC_OUT}", ' '.join([cor(C.__isoform_dense_star_splice_output__)] + cor(C.__isoform_dense_star_splice_output_logs__)), 'pipeline_isoform_dense_star_splice.py');
    M.add_step("ISOFORM_DENSE_STAR_GG", "${ISOFORM_DENSE_SPLIT_GENOME_OUT} ${ISOFORM_DENSE_STAR_SPLICE_OUT}", C.__isoform_dense_genome_generate_dir__(), 'pipeline_isoform_dense_star_genome_generate.py');
  else:
    M.add_step("ISOFORM_DENSE_STAR_GG", "${ISOFORM_DENSE_SPLIT_GENOME_OUT}", C.__isoform_dense_genome_generate_dir__(), 'pipeline_isoform_dense_star_genome_generate.py');
  #fi
  M.add_var("ISOFORM_DENSE_STAR_UNMAPPED", ' '.join(flatten(C.__isoform_dense_star_align_output_unmapped__())));
  M.add_var("ISOFORM_DENSE_STAR_BAM", ' '.join(C.__isoform_dense_star_align_output_bam__()));
  M.add_var("ISOFORM_DENSE_STAR_LOGS", ' '.join(flatten(C.__isoform_dense_star_align_output_log__())));
  M.add_var("ISOFORM_DENSE_STAR_MERGED_BAM", C.__isoform_dense_star_align_output_merged__());
  M.add_step("ISOFORM_DENSE_STAR", "${ISOFORM_DENSE_STAR_GG_OUT} ${TRIMMOMATIC_OUT}", '${ISOFORM_DENSE_STAR_UNMAPPED} ${ISOFORM_DENSE_STAR_BAM} ${ISOFORM_DENSE_STAR_LOGS} ${ISOFORM_DENSE_STAR_MERGED_BAM} ', 'pipeline_isoform_dense_star_align.py');
  M.add_step("ISOFORM_DENSE_CUFFLINKS_RABT", "${ISOFORM_DENSE_STAR_MERGED_BAM} ${ISOFORM_DENSE_SPLIT_GENOME_OUT}", ' '.join(cor(C.__isoform_dense_cufflinks_output__)),  'pipeline_isoform_dense_cufflinks_rabt.py');
  M.add_step("ISOFORM_DENSE_UNSPLIT_GENOME", "${ISOFORM_DENSE_CUFFLINKS_RABT_OUT}", cor(C.__isoform_dense_genome_unsplit_output__), 'pipeline_isoform_dense_genome_unsplit.py');
  M.add_step("ISOFORM_DENSE_UNSPLIT_ANALYSIS", "${ISOFORM_DENSE_UNSPLIT_GENOME_OUT}", cor(C.__isoform_dense_genome_analysis_outdir__), 'pipeline_isoform_dense_analysis.py');
#fi

M.write(C.makefile, C.location, "${QUALITYREPORT} ${POSTANALYSIS}", C.outdir, C.jobname);

###############################################################################

sys.exit(0);
