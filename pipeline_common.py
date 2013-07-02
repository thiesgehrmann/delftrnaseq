#!/usr/bin/python

import os;
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

  email          = None;
  jobname        = None;
  workdir        = "./"
  outdir         = "./"
  makefile       = "Makefile";
  location       = "";
  samples        = [];
  sample_names   = [];
  sample_labels  = [];
  sample_comp    = [];
  genome         = [];
  genome_annot   = None;
  genome_guide   = None;

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

  trimmomatic_opts="-threads 12";
  trimmomatic_trim="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";

  def __trimmomatic_output__(self):
    return [ (self.outdir + '/' + self.sample_names[i] + '_R1.fastq', self.outdir + '/' + self.sample_names[i] + '_R2.fastq') for i in xrange(len(self.sample_names)) ];
  #edef

  #############################################################################
  # STAR GENOME GENERATION STUFF                                              #
  #############################################################################

  def star_gg_opts(self):
    return "--runThreadN %d" % self.__max_threads__;
  #edef

  star_gg_gdir = "genome"

  def __star_gg_output__(self):
    return "%s/%s" %(self.outdir, self.star_gg_gdir);
  #edef

  #############################################################################
  # STAR ALIGNMENT STUFF                                                      #
  #############################################################################

  def star_al_opts(self):
    return "--runThreadN %d" % self.__max_threads__;
  #edef

  def __star_al_output_sam__(self):
    return [ self.outdir + '/%s.star_align.sam' % sn for sn in self.sample_names ];
  #edef
  def __star_al_output_unmapped__(self):
    return [ (self.outdir + "/%s.star_align_unmapped_R1.fastq" % sn, self.outdir + "/%s.star_align_unmapped_R2.fastq" % sn) for sn in self.sample_names ];
  #edef

  #############################################################################
  # POST STAR AL BAM STUFF                                                    #
  #############################################################################

  def __post_star_al_bam_output__(self):
    return [ self.outdir + "/%s.star_align.bam" % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # POST STAR AL SORT STUFF                                                   #
  #############################################################################

  def __post_star_al_sort_output__(self):
    return [ self.outdir + "/%s.star_align_sort.bam" % sn for sn in self.sample_names ];
  #edef

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

  def cuffdiff_opts(self):
    return "-p %s --upper-quartile-norm" % self.__max_threads__;
  #edef

  def __cuffdiff_output__(self):
    files = [ "bias_params.info", "cds.count_tracking", "cds.diff", "cds.fpkm_tracking", "cds.read_group_tracking", "cds_exp.diff", "gene_exp.diff", "genes.count_tracking", "genes.fpkm_tracking", "genes.read_group_tracking", "isoform_exp.diff", "isoforms.count_tracking", "isoforms.fpkm_tracking", "isoforms.read_group_tracking", "promoters.diff", "read_groups.info", "run.info", "splicing.diff", "tss_group_exp.diff", "tss_groups.count_tracking", "tss_groups.fpkm_tracking", "tss_groups.read_group_tracking", "var_model.info" ];
    r = [];
    for (a,b) in self.cuffdiff_cmp:
      r.append(tuple([ "%s/%s-%s,%s.cuffdiff.%s" % (self.outdir, self.jobname, self.label_names[a], self.label_names[b], f) for f in files ]));
    #efor
    return r;
  #edef

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
  def unmapped_blast_opts():
    return "-num_threads %d" % __max_threads__;
  #edef

  __unmapped_blast_fields__  = "qseqid sseqid slen length mismatch gapopen pident evalue bitscore";

  def __unmapped_blast_output__(self):
    return [ self.outdir + '/%s.unmapped_blast.tsv' % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # UNMAPPED STUFF                                                            #
  #############################################################################
  
  unmapped_opts = "";
  __unmapped_blast_select_by__ = len(__unmapped_blast_fields__.split(' ')) - 1;

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

    self.genome = wd + '/' + genome;

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

      # Check if sample files are there
    for (r1, r2) in self.samples:
      errors = errors + fex(r1, "Could not find sample file '%s'" % r1);
      errors = errors + fex(r2, "Could not find sample file '%s'" % r2);
    #efor

    if warnings > 0:
      print("There were %d warnings!" % warnings);
    #fi
    if errors > 0:
      print("There were %d errors! Aborting!" % errors);
      return 1
    #fi
    return 0
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

def run_cmd(cmd, bg=False, stdin=None, stdout=None, stderr=None):
  p = subprocess.Popen(shlex.split(cmd), stdin=stdin, stdout=stdout, stderr=stderr);
  if bg:
    return p;
  else:
    (pid, r) = os.waitpid(p.pid, 0);
    return r;
  #fi
#edef

###############################################################################

def run_par_cmds(cmd_list, max_threads=PIPELINECONF.__max_threads__, stdin=None, stdout=None, stderr=None):
  
  p = [];
  i = 0;
  retval = 0;
  cmds = len(cmd_list);

  while i < cmds:
    while len(p) < max_threads and i < cmds:
      print "RUNNING: %s" % cmd_list[i]; sys.stdout.flush();
      p.append( (run_cmd(cmd_list[i], bg=True, stdin=stdin, stdout=stdout, stderr=stderr),i) );
      i = i + 1;
    #ewhile

    time.sleep(0.5);

    running   = [ (j, k) for (j,k) in p if j.poll() == None ];
    completed = [ (j, k) for (j,k) in p if j.poll() != None ];

    for (j,k) in completed:
      if j.returncode != 0:
        retval = retval + j.returncode;
        print "ERROR: Failed in cmd: %s" % cmd_list[k]; sys.stdout.flush();
      else:
        print "COMPLETED: cmd : %s" % cmd_list[k]; sys.stdout.flush();
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

