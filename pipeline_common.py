#!/usr/bin/python

import os
import types
import pickle
import tempfile
import shlex, subprocess
###############################################################################

class PIPELINECONF:

  inst_loc = None;

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
  genome         = [];
  genome_annot   = None;
  genome_guide   = None;

  #############################################################################
  # INTERNALS                                                                 #
  #############################################################################

  __pipeline_email__       = 'delftrnaseq@gmail.com';
  __pipeline_mail_user__   = 'delftrnaseq';
  __pipeline_mail_pass__   = 'thiesgehrmann';
  __pipeline_mail_server__ = 'smtp.gmail.com:587';
 


  #############################################################################
  # TRIMM-O-MATIC STUFF                                                       #
  #############################################################################

  trimmomatic_opts="-threads 12";
  trimmomatic_trim="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";

  def trimmomatic_output(self):
    return [ (self.outdir + '/' + self.sample_names[i] + '_R1.fastq', self.outdir + '/' + self.sample_names[i] + '_R2.fastq') for i in xrange(len(self.sample_names)) ];
  #edef

  #############################################################################
  # STAR GENOME GENERATION STUFF                                              #
  #############################################################################

  star_gg_opts = "--runThreadN 16";
  star_gg_gdir = "genome"

  def star_gg_output(self):
    return "%s/%s" %(self.outdir, self.star_gg_gdir);
  #edef

  #############################################################################
  # STAR ALIGNMENT STUFF                                                      #
  #############################################################################

  star_al_opts="--runThreadN 16";

  def star_al_output_sam(self):
    return [ self.outdir + '/%s.star_align.sam' % sn for sn in self.sample_names ];
  #edef
  def star_al_output_unmapped(self):
    return [ (self.outdir + "/%s.star_align_unmapped_R1.fastq" % sn, self.outdir + "/%s.star_align_unmapped_R2.fastq" % sn) for sn in self.sample_names ];
  #edef

  #############################################################################
  # POST STAR AL STUFF                                                        #
  #############################################################################

  def post_star_al_output(self):
    return [ self.outdir + "/%s.star_align_sort_name.bam" % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # GENOME GENERATION STUFF                                                   #
  #############################################################################

  def genome_gen_output(self):
    return self.outdir + '/genome.gff';
  #edef

  #############################################################################
  # PRE-CUFFLINKS STUFF                                                       #
  #############################################################################
  
  pre_cufflinks_opts="";
  
  def pre_cufflinks_output(self):
    return self.outdir + '/cufflinks_prep.bam';
  #edef

  #############################################################################
  # CUFFLINKS STUFF                                                           #
  #############################################################################

  cufflinks_opts="-p 12 -u --max-intron-length 5000 --min-intron-length 25 --overlap-radius 25 --max-bundle-length 250000";
  cufflinks_bias_corr="";

  def cufflinks_output(self):
    return [ self.outdir + f for f in [ '/transcripts.gtf', '/skipped.gtf', '/isoforms.fpkm_tracking', '/genes.fpkm_tracking' ] ]
  #edef

  #############################################################################
  # CUFFDIFF STUFF                                                            #
  #############################################################################

  cuffdiff_opts = "";
  def cuffdiff_output(self):
    return ["somethingidontknow_cd"];
  #edef

  #############################################################################
  # TRINITY STUFF                                                             #
  #############################################################################

  trinity_opts = "--CPU 6 --jaccard_clip --min_contig_length=100 -bflyCalculateCPU";

  def trinity_output(self):
    return [ self.outdir + "/%s_trinity_assembled.fasta" % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # TRINITY_ORF STUFF                                                         #
  #############################################################################

  trinity_orf_opts = "";
  def trinity_orf_output(self):
    return [ self.outdir + '/%s.trinity_orfs.fasta' % sn for sn in self.sample_names ];
  #edef

  #############################################################################
  # UNMAPPED STUFF                                                            #
  #############################################################################
  
  unmapped_opts = "";
  def unmapped_output(self):
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

    genome = []
    for f in self.genome:
      genome = genome + [ wd + '/' + f ];
    #efor
    self.genome = genome

    self.genome_annot = wd + '/' + self.genome_annot if self.genome_annot else None;
    self.genome_guide = wd + '/' + self.genome_guide if self.genome_guide else None;
  #edef

  #############################################################################

  def data_check(self):

    errors   = 0
    warnings = 0

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
    print p;


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

def run_par_cmds(cmd_list):

  for cmd in [ x for x in cmd_list if x ]:
    subprocess.Popen(shlex.split(cmd));
  #efor
  for cmd in [ x for x in cmd_list if x ]:
    os.wait();
  #efor
#edef

###############################################################################

def run_cmd(cmd, stdin=None, stdout=None, stderr=None):
  return subprocess.call(shlex.split(cmd), stdin=stdin, stdout=stdout, stderr=stderr)
#edef

###############################################################################

def run_seq_cmds(cmd_list):
  for cmd in [ x for x in cmd_list if x ]:
    retval = run_cmd(cmd);
    if retval != 0:
      return retval;
    #fi
  #efor

  return 0;

###############################################################################

def flatten(l):
  import collections

  for el in l:
    if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
      for sub in flatten(el):
        yield sub
      #efor
    else:
      yield el
    #fi
  #efor
#edef

###############################################################################

