#!/usr/bin/env python

import os;
import pwd;
import sys;
import types;
import pickle;
import tempfile;
import shlex, subprocess;
import multiprocessing;
import time;

###############################################################################

class PIPELINECONF:

  #############################################################################
  # DATA                                                                      #
  #############################################################################

  email           = None;
  jobname         = None;
  title           = "Rna-Seq Analysis";
  author          = pwd.getpwuid(os.getuid()).pw_gecos;
  
  workdir         = "./"
  outdir          = "./"
  makefile        = "Makefile";

  location        = "";
  PE              = None;
  strand_specific = False;
  samples         = [];
  sample_names    = [];
  sample_labels   = [];
  sample_comp     = [];
  genome          = [];
  genome_annot    = None;
  genome_guide    = None;

  #############################################################################
  # INTERNALS                                                                 #
  #############################################################################

  inst_loc    = None;
  __max_threads__ = multiprocessing.cpu_count();
  __pipeline_email__       = 'delftrnaseq@gmail.com';
  __pipeline_mail_user__   = 'delftrnaseq';
  __pipeline_mail_pass__   = 'thiesgehrmann';
  __pipeline_mail_server__ = 'smtp.gmail.com:587';
 


  #############################################################################
  # TRIMM-O-MATIC STUFF                                                       #
  #############################################################################

  def trimmomatic_opts(self):
    return "-threads %d" % self.__max_threads__;
  #edef

  trimmomatic_trim="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";

  def __trimmomatic_output__(self):
    if self.PE:
      return [ (self.outdir + '/' + sn + '_R1.fastq', self.outdir + '/' + sn + '_R2.fastq') for sn in self.sample_names ];
    else:
      return [ self.outdir + '/' + sn + '_READS.fastq' for sn in self.sample_names ];
    #fi
  #edef

  #############################################################################
  # STAR GENOME GENERATION STUFF                                              #
  #############################################################################
  def star_pregg_opts(self):
    return "--runThreadN %d" % self.__max_threads__;
  #edef

  def star_gg_opts(self):
    return "--runThreadN %d --sjdbOverhang 99" % self.__max_threads__;
  #edef

    # Do we build the spliceDB?
  build_splice_db = True;

  star_pre_splice_gdir = "pre_splice_genome"
  def __star_pre_splice_output__(self):
    return "%s/%s" %(self.outdir, self.star_pre_splice_gdir);
  #edef

  star_gg_gdir = "genome"
  def __star_gg_output__(self):
    return "%s/%s" %(self.outdir, self.star_gg_gdir);
  #edef

  #############################################################################
  # STAR ALIGNMENT STUFF                                                      #
  #############################################################################
  def star_preal_opts(self):
    return "--runThreadN %d" % self.__max_threads__;

  #splice sites filtering
  #condition: (min_overhang >= $SPLICE_SITES_MINMAX_OVERHANG and at least 2 unique reads) OR unique_reads > $SPLICE_SITES_UNIQUE_READS OR nonunique_reads > SPLICE_SITES_NONUNIQUE_READS
  splice_sites_minmax_overhang = 36 #(and minimal 2 unique reads) #OR
  def splice_sites_unique_reads(self): #OR
    return round(len(self.samples) / 2)
  def splice_sites_nonunique_reads(self): #OR
    return len(self.samples)

  def star_al_opts(self):
    return "--runThreadN %d" % self.__max_threads__;
  #edef
  
  def __star_preal_output__(self):
    return self.outdir + "/splice_junction_db_nohead.tsv"

  def __star_al_output_sam__(self):
    return [ self.outdir + '/%s.star_align.sam' % sn for sn in self.sample_names ];
  #edef
  def __star_al_output_unmapped__(self):
    if self.PE:
      return [ (self.outdir + "/%s.star_align_unmapped_R1.fastq" % sn, self.outdir + "/%s.star_align_unmapped_R2.fastq" % sn) for sn in self.sample_names ];
    else:
      return [ self.outdir + "/%s.star_align_unmapped_R1.fastq" % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # POST STAR AL BAM STUFF                                                    #
  #############################################################################

  max_bam_threads = min(__max_threads__ / 2, 4) #reduce number of threads for bam (are a bit intensive on the disk)


  def __post_star_al_bam_output__(self) :
    return [ self.outdir + "/%s.star_align.bam" % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # POST STAR AL SORT STUFF                                                   #
  #############################################################################

  def __post_star_al_sort_output__(self):
    return [ self.outdir + "/%s.star_align_sort.bam" % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # POST STAR AL INDEX STUFF                                                  #
  #############################################################################

  def __post_star_al_index_output__(self):
    return [ self.outdir + "/%s.star_align_sort.bam.bai" % sn for sn in self.sample_names ];
  
  #############################################################################
  # CDS GFF STUFF                                                             #
  #############################################################################

  cds_gff_type = 'CDS'

  def __cds_gff_output__(self):
    return [ self.outdir + "/%s.star_align_sort.count" % sn for sn in self.sample_names ]

  #############################################################################
  # READ DISTRIBUTION STUFF                                                   #
  #############################################################################

  def __read_distribution_output__(self):
    return [ self.outdir + "/%s.star_align_sort.read_stats.pdf" % sn for sn in self.sample_names ]

  #############################################################################
  # CDS GFF STUFF                                                             #
  #############################################################################
  
  def __fastqc_output__(self) : 
    out = []
    for bam_output in self.__post_star_al_bam_output__() :
        prefix, ext = os.path.splitext(bam_output)
        filename = prefix + '_fastqc'
        out.append(filename)
    if self.PE :
        for left_in, right_in in self.__star_al_output_unmapped__() :
            prefix, ext = os.path.splitext(left_in)
            out.append(prefix + '_fastqc')
            prefix, ext = os.path.splitext(right_in)
            out.append(prefix + '_fastqc')
    else :
        for file in self.__star_al_output_unmapped__() :
            prefix, ext = os.path.splitext(file)
            output.append(prefix + '_fastqc')
            out.append(prefix + '_fastqc')
    return out

  #############################################################################
  # GENOME GENERATION STUFF                                                   #
  #############################################################################

  def __genome_gen_output__(self):
    return self.outdir + '/genome.gff';
  #edef

  #############################################################################
  # GENOME ANNOT FORMAT STUFF                                                 #
  #############################################################################

  def __genome_annot_format_output__(self):
    return "%s/%s.cleaned.gff" % (self.outdir, self.jobname);
  #edef

  #############################################################################
  # PRE-CUFFLINKS MERGE STUFF                                                 #
  #############################################################################
  
  pre_cufflinks_merge_opts="";
  
  def __pre_cufflinks_merge_output__(self):
    return self.outdir + '/%s.pre_cufflinks_merged.bam' % self.jobname;
  #edef

  #############################################################################
  # PRE-CUFFLINKS SORT STUFF                                                  #
  #############################################################################

  pre_cufflinks_sort_opts="";

  def __pre_cufflinks_sort_output__(self):
    return self.outdir + '/%s.pre_cufflinks_sorted.bam' % self.jobname;
  #edef

  #############################################################################
  # CUFFLINKS STUFF                                                           #
  #############################################################################

  def cufflinks_opts(self):
    return "-p %d -u --max-intron-length 5000 --min-intron-length 25 --overlap-radius 25 --max-bundle-length 250000" % self.__max_threads__;
  #edef

  cufflinks_bias_corr="";

  def __cufflinks_output__(self):
    return [ "%s/%s.cufflinks.%s" % (self.outdir, self.jobname, f) for f in [ "genes.fpkm_tracking", "isoforms.fpkm_tracking", "transcripts.gtf", "skipped.gtf" ] ];
  #edef

  #############################################################################
  # PRE-CUFFLINKS_INDIV SORT STUFF                                            #
  #############################################################################

  pre_cufflinks_indiv_sort_opts="";

  def __pre_cufflinks_indiv_sort_output__(self):
    return [ self.outdir + '/%s.pre_cufflinks_indiv_sorted.bam' % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # CUFFLINKS INDIV STUFF                                                     #
  #############################################################################

  def cufflinks_indiv_opts(self):
    return "-p %d -u --max-intron-length 5000 --min-intron-length 25 --overlap-radius 25 --max-bundle-length 250000" % self.__max_threads__;
  #edef

  cufflinks_bias_corr="";

  def __cufflinks_indiv_output__(self):
    r = [];
    for sn in self.sample_names:
      r.append([ "%s/%s.cufflinks_indiv.%s" % (self.outdir, sn, f) for f in [ "genes.fpkm_tracking", "isoforms.fpkm_tracking", "transcripts.gtf", "skipped.gtf" ] ]);
    #efor
    return r;
  #edef

  #############################################################################
  # CUFFDIFF STUFF                                                            #
  #############################################################################

  cuffdiff_cmp = None;
  cuffdiff_test_type         = 'isoform'; # cds|gene|isoform|tss
  __cuffdiff_test_group_index__  = {'cds':3, 'gene':9, 'isoform':13, 'tss':21};
  __cuffdiff_test_diff_index__   = {'cds':5, 'gene':6, 'isoform':10, 'tss':18};

  def cuffdiff_opts(self):
    return "-p %s --upper-quartile-norm --max-bundle-frags 100000000000" % self.__max_threads__;
  #edef

  def __cuffdiff_output__(self):
    files = [ "bias_params.info", "cds.count_tracking", "cds.diff", "cds.fpkm_tracking", "cds.read_group_tracking", "cds_exp.diff", "gene_exp.diff", "genes.count_tracking", "genes.fpkm_tracking", "genes.read_group_tracking", "isoform_exp.diff", "isoforms.count_tracking", "isoforms.fpkm_tracking", "isoforms.read_group_tracking", "promoters.diff", "read_groups.info", "run.info", "splicing.diff", "tss_group_exp.diff", "tss_groups.count_tracking", "tss_groups.fpkm_tracking", "tss_groups.read_group_tracking", "var_model.info" ];
    return [ "%s/%s-all.cuffdiff.%s" % (self.outdir, self.jobname, f) for f in files ];
  #edef
  
  #############################################################################
  # CUFFDIFF_COMBINE STUFF                                                    #
  #############################################################################

  annotation_files = [];
  annotation_names = [];
  def __cuffdiff_combine_output__(self):
    return [ self.outdir + '/cuffdiff_combine.dat', self.outdir + '/cuffdiff_combine.csv'];
  #edef

  #############################################################################
  # ANALYSIS STUFF                                                            #
  #############################################################################

  perform_quality_report     = True;
  perform_analysis           = True;
  analysis_filter            = [None];
  analysis_filter_names      = [ "NoFilter" ];
  analysis_venn_updown_split = False;
  analysis_venn              = [];

  analysis_enrichment_verbose_output  = False;
  analysis_enrichment_returned_values = None;
  analysis_enrichment_only_annotated  = False; # Count only the annotated genes, or all genes (set to true if you set updown_split to true)
  analysis_enrichment_updown_split    = False; # Enrich for sets of differentially up and diferentially downregulated genes.
  __analysis_enrichment_alpha__       = 0.05; 

  def __analysis_output__(self):
    files = []
    for filter_name in self.analysis_filter_names:
      filtfiles = [];
      if self.analysis_venn_updown_split:
        filtfiles.append([['%s/%s_filter=%s_venn_%d_up.pdf' % (self.outdir, self.jobname, filter_name, pos), '%s/%s_filter=%s_venn_%d_down.pdf' % (self.outdir, self.jobname, filter_name, pos) ]for pos in range(len(self.analysis_venn))])
      else:
        filtfiles.append([['%s/%s_filter=%s_venn_%d.pdf' % (self.outdir, self.jobname, filter_name, pos)] for pos in range(len(self.analysis_venn))])
      #fi
      filtfiles.append(['%s/%s_filter=%s_pca.pdf' % (self.outdir, self.jobname, filter_name), '%s/%s_filter=%s_pca_log_filter.pdf' % (self.outdir, self.jobname, filter_name)])
      filtfiles.append(['%s/%s_filter=%s_diffstats.pdf' % (self.outdir, self.jobname, filter_name)])
      filtfiles.append(['%s/%s_filter=%s_clusters.pdf' % (self.outdir, self.jobname, filter_name)])
      filtfiles.append(['%s/%s_filter=%s_venn_subsets.csv' % (self.outdir, self.jobname, filter_name), '%s/%s_filter=%s_venn_subsets.dat' % (self.outdir, self.jobname, filter_name)]);
      filtfiles.append(['%s/%s_filter=%s_enrichment.dat' % (self.outdir, self.jobname, filter_name)]);
      filtfiles.append(['%s/%s_filter=%s_analysis_report.tex' % (self.outdir, self.jobname, filter_name)])
      files.append(filtfiles);
    #efor
    return files;
  #edef


  def __quality_output__(self):
    files = [];
    if self.PE:
      files.append(['%s/%s_trimmomatic_nreads.pdf' % (self.outdir, self.jobname), '%s/%s_trimmomatic_ratio.pdf' % (self.outdir, self.jobname), '%s/%s_trimmomatic_single.pdf' % (self.outdir, self.jobname)])
    else:
      files.append(['%s/%s_trimmomatic_nreads.pdf' % (self.outdir, self.jobname), '%s/%s_trimmomatic_ratio.pdf' % (self.outdir, self.jobname)]);
    #fi
    files.append(['%s/%s_mapping_matched.pdf' % (self.outdir, self.jobname), '%s/%s_mapping_unmatched.pdf' % (self.outdir, self.jobname), '%s/%s_mapping_length.pdf' % (self.outdir, self.jobname)])
    files.append(['%s/%s_quality_report.tex' % (self.outdir, self.jobname)])
    return files
  #edef


  #############################################################################
  # CONTAMINATION STUFF                                                       #
  #############################################################################

  check_contamination = False;

  #############################################################################
  # TRINITY STUFF                                                             #
  #############################################################################

  def trinity_opts(self):
    return "--CPU %d --jaccard_clip --min_contig_length=100 -bflyCalculateCPU" % self.__max_threads__;
  #edef

  def __trinity_output__(self):
    return [ self.outdir + "/%s.trinity_assembled.fasta" % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # TRINITY_ORF STUFF                                                         #
  #############################################################################

  trinity_orf_opts = "";

  def __trinity_orf_output__(self):
    return [ self.outdir + '/%s.trinity_orfs.fasta' % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # UNMAPPED BLAST STUFF                                                      #
  #############################################################################

  blast_db = None;
  def unmapped_blast_opts(self):
    return "-num_threads %d" % self.__max_threads__;
  #edef

  __unmapped_blast_fields__  = "qseqid sseqid qlen slen length mismatch gapopen pident evalue bitscore";

  def __unmapped_blast_output__(self):
    return [ self.outdir + '/%s.unmapped_blast.tsv' % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # UNMAPPED STUFF                                                            #
  #############################################################################
  
  unmapped_opts = "";
  
  __unmapped_blast_select_by__ = len(__unmapped_blast_fields__.split(' ')) - 1;
  
  __unmapped_blast_select_by_cutoff__ = 750
  
  unmapped_use_mysql = False
  unmapped_mysql_user = 'genbank'
  unmapped_mysql_pass = 'genbank'
  unmapped_mysql_host = 'localhost'
  unmapped_mysql_db = 'genbank'

  def __unmapped_output__(self):
    return [ self.outdir + '/%s.unmapped_orgs.dat' % sn for sn in self.sample_names ];
  #edef
  
  #############################################################################
  #############################################################################  
  #############################################################################
  #############################################################################
  #############################################################################
  #############################################################################

  def __init__(self, pname = None):
    if not(pname is None):
      # Unpickle it!
      p = pickle.load(open(pname, 'r'));
      for (attr, val) in p:
        setattr(self, attr, val);
      #efor
    #fi
  #edef

  #############################################################################

  def wd_correct(self):
    wd = self.workdir

    if wd == None:
      return
    #fi

    samples = []
    for (r1, r2) in self.samples:
      samples = samples + [ (wd + '/' + r1, wd + '/' + r2) ];
    #efor
    self.samples = samples

    self.genome = wd + '/' + self.genome;

    self.genome_annot = wd + '/' + self.genome_annot if self.genome_annot else None;
    self.genome_guide = wd + '/' + self.genome_guide if self.genome_guide else None;
  #edef

  #############################################################################

  def data_check(self):

    errors   = 0
    warnings = 0

    if self.jobname == None or self.jobname == "":
      self.jobname = 'UNDEFINED';
      warnings = warnings + 1;
      warning("No job name given, please set 'jobname' appropriately. Using '%s'." % self.jobname);
    #fi

    if (self.genome_annot == None) and (self.genome_guide == None):
      warnings = warnings + 1;
      warning("No Genome annotation given, please set 'genome_annot' or 'genome_guide' variables");
    #fi

    if self.genome_annot:
      errors = errors + fex(self.genome_annot, "Could not find genome annotation '%s'" % self.genome_annot);
    #fi
    if self.genome_guide:
      errors = errors + fex(self.genome_guide, "Could not find genome guide '%s'" % self.genome_guide);

    if not len(self.samples) == len(self.sample_names):
      warnings = warnings + 1
      self.sample_names = [ "%d" %(i+1) for i in xrange(len(self.samples)) ]
      warning("Samples names are set incorrectly. Please set 'sample_names' variable.");
    #fi

    if not len(self.sample_labels) == len(self.samples):
      warnings = warnings + 1;
      self.sample_labels = self.sample_names;
      warning("Sample labels are set incorrectly. Please set 'sample_labels' variable.");

    if len(self.samples) == 0:
      error("No samples specified; What do you expect me to do??");
      errors = errors + 1;
    else:
      if self.PE == None:
        self.PE = hasattr(self.samples[0], '__iter__');
        warning("Please specify if you are using paired end or single end reads, set variable PE. Detecting %s." % str(self.PE));
        warnings = warnings + 1;
        # Check if sample files are there
      #fi
      for rs in self.samples:
        if self.PE:
          errors = errors + fex(rs[0], "Could not find sample file '%s'" % rs[0]);
          errors = errors + fex(rs[1], "Could not find sample file '%s'" % rs[1]);
        else:
          errors = errors + fex(rs, "Could not find sample file '%s'" % rs);
        #fi
      #efor
    #fi

    n_cmps = len(self.cuffdiff_cmp);
    self.cuffdiff_cmp = list(set(self.cuffdiff_cmp));
    if len(self.cuffdiff_cmp) < n_cmps:
      warnings = warnings + 1;
      warning("Some cuffdiff comparisons are specified more than once, correcting. May cause confusion in later analysis");
    #fi

    for file in self.annotation_files:
      errors = errors + fex(file, "Could not find annotation file '%s'." % file);
    #efor
    if not(len(self.annotation_files) == len(self.annotation_names)):
      errors = errors + 1;
      error("The annotation files in the 'annoation_files' variable have not been given names properly. Please set 'annotation_names' variable.");
    #fi

    if self.perform_analysis:

      if not(len(self.analysis_filter)  == len(self.analysis_filter_names)):
        warnings = warnings + 1;
        warning("Filters have not been assigned names properly. Please set analysis_filter_names correctly. Assigning names.");
        self.analysis_filter_names = [ str(f) for f in self.analysis_filter ];
      #fi

      for filt in self.analysis_filter:
        if filt == None:
          continue;
        #fi
        errors = errors + fex(filt, "Could not find filter file '%s'." % filt);
      #efor

      for venn in self.analysis_venn:
        if len(venn) not in [2,3]:
          errors = errors + 1;
          error("The venn analysis: '%s' is impossible. Need 2 or 3 cuffdiff comparisons." % (str(venn)));
        #fi
        for v in venn:
          if v < 0 or v > len(self.cuffdiff_cmp):
            errors = errors + 1;
            error("The venn analysis: '%s' is impossible. Cuffdiff comparison %d does not exist." % (str(venn), v));
          #fi
        #efor
      #efor

    #fi

    if warnings > 0:
      print("There were %d warnings!" % warnings);
    #fi
    if errors > 0:
      print("There were %d errors! Aborting!" % errors);
      return 1;
    #fi
    return 0;
  #edef

  #############################################################################

  def pickle(self):
    isfunc = lambda obj, attr: hasattr(obj, attr) and type(getattr(obj, attr)) == types.MethodType;
    p = [ (attr, getattr(self, attr)) for attr in dir(self) if not(isfunc(self, attr)) ];

    fd = None;
    if self.location == "":
      fd, fname = tempfile.mkstemp();
      p = p + [ ('location', fname) ];
    else:
      fname = self.location;
    #fi
    f = open(fname, 'w');
    pickle.dump(p, f);
    f.close();
    if fd:
      os.close(fd);
    #fi
    return fname;
  #edef

#eclass



def usage(a1):
    print "Usage:  %s <config file>" % a1;

def init_conf():
    if len(os.sys.argv) != 2:
        usage(os.sys.argv[0]);
        os.sys.exit(1);
    
    C = PIPELINECONF(os.sys.argv[1])
    run_cmd('mkdir -p %s' % C.outdir);
    return C
      
   


###############################################################################
###############################################################################
###############################################################################
###############################################################################

def write_config(c):
  import pickle
  import tempfile

  if c.conffile == "":
    fd, temp_path = tempfile.mkstemp();
    f = open(temp_path, 'w');
  else:
    print c.conffile
    fd = os.fdopen(c.conffile, 'w');
    f = open(c.conffile, 'w');
  pickle.dump(c, f);
  f.close();
  os.close(fd);

#edef

###############################################################################

def fex(filename, str):
  if os.path.isfile(filename):
    return 0
  #fi
  error(str)
  return 1
#edef

###############################################################################

def fnex(filename, str):
  if not(filename) or os.path.isfile(filename):
    return 0
  #fi
  warning(str);
  return 1;
#edef

###############################################################################

def warning(str):
  os.sys.stderr.write("[W]" + str + '\n');
#edef

###############################################################################

def error(str):
  os.sys.stderr.write("[E]" + str + '\n');
#edef

###############################################################################

def cor(obj):
  """ Call or return.
      Returns obj() if obj is callable, otherwise returns obj.
  """

  if hasattr(obj, '__call__'):
    return obj();
  else:
    return obj;
#edef

###############################################################################

def run_cmd(cmd, bg = False, stdin = None, stdout = None, stderr = None):
    print '[RUN CMD] Executing: %s' % cmd
    sys.stdout.flush()
    p = subprocess.Popen(shlex.split(cmd), stdin=stdin, stdout=stdout, stderr=stderr);
    if bg :
        return p
    else:
        (pid, r) = os.waitpid(p.pid, 0);
        return r;

def run_shell(cmd, bg = False):
    print '[RUN SHELL] Executing: %s' % cmd
    sys.stdout.flush()
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    if bg :
        return p
    else:
        (pid, r) = os.waitpid(p.pid, 0);
        return r;

def getCommandOutput(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    if p.returncode != 0:
        raise RuntimeError, "Command failed: " + str(res)
    return res

###############################################################################

def run_par_cmds(cmd_list, max_threads=PIPELINECONF.__max_threads__, stdin=None, stdout=None, stderr=None, shell=False):
  
  p = [];
  i = 0;
  retval = 0;
  cmds = len(cmd_list);

  while i < cmds:
    while len(p) < max_threads and i < cmds:
      print "RUNNING: %s" % cmd_list[i]; sys.stdout.flush();
      sys.stdout.flush()
      if not shell :
        p.append( (run_cmd(cmd_list[i], bg=True, stdin=stdin, stdout=stdout, stderr=stderr),i) );
      else :
        p.append( (run_shell(cmd_list[i], bg=True),i) );
      i = i + 1;
    #ewhile

    time.sleep(0.5);

    running   = [ (j, k) for (j,k) in p if j.poll() == None ];
    completed = [ (j, k) for (j,k) in p if j.poll() != None ];

    for (j,k) in completed:
      if j.returncode != 0:
        retval = retval + j.returncode;
        print "ERROR: Failed in cmd: %s" % cmd_list[k]; sys.stdout.flush();
        sys.stdout.flush()
      else:
        print "COMPLETED: cmd : %s" % cmd_list[k]; sys.stdout.flush();
        sys.stdout.flush()
      #fi
    #efor
    p = running;
  #ewhile
  
  return retval;
#edef

###############################################################################

def run_seq_cmds(cmd_list, stdin=None, stdout=None, stderr=None):

  for cmd in [ x for x in cmd_list if x ]:
    retval = run_cmd(cmd, stdin=stdin, stdout=stdout, stderr=stderr);
    if retval != 0:
      print "ERROR: Failed on cmd: %s" % cmd;
      sys.stdout.flush()
      return retval;
    #fi
  #efor

  return 0;
#edef

###############################################################################

def flatten(l):
  import collections

  for el in l:
    if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
      for sub in flatten(el):
        yield sub;
      #efor
    else:
      yield el;
    #fi
  #efor
#edef

###############################################################################

