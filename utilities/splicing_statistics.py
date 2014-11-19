import sys;
from bx.intervals.intersection import IntervalTree
import numpy as np;
from ibidas import *;
import matplotlib.pylab as plt;

###############################################################################

def index_gff3_id(G, features=["exon"]):

  G = G[_.feature.In(Rep(features))];
  #G = G[_.strand == '+'].GroupBy(_.parent).Sort(_.start, descend=False) | Stack | G[_.strand == '-'].GroupBy(_.parent).Sort(_.start, descend=True);
  G = G.GroupBy(_.parent).Sort(_.start);
  G = G.Get(_.parent, _.start, _.end, _.strand).Flat();

  genes = {};
  for (parent, start, end, strand) in zip(*G()):
    if parent not in genes:
      genes[parent] = [{}, IntervalTree()];
    #fi

    exon_n = len(genes[parent][0].keys());
    genes[parent][0][exon_n+1] = (start, end);
    genes[parent][1].add(start, end, (start, end, exon_n+1));
    print parent, start, end, exon_n+1;
  #efor
    
  return genes;
#edef

###############################################################################
# Diagnose a given alternative splicing construct

# 0: Region Start
# 1: Region End
# 2: Mapped exons
# 3: Intron retention
# 4: Alternative 5' SS
# 5: Alternative 3' SS
# 6: Exon skipping
# 7: New exon
# 8: Identical exon
diag_indexes = [ 0,1,2,3,4,5,6,7,8,9,10 ];
RS, RE, ORS, ORE, MP, IR, A5, A3, ES, NE, IE = diag_indexes;

def diagnose_exon(exon_range, gene_exon_index):

  diagnosis = [False]*len(diag_indexes);
  diagnosis[RS] = exon_range[0];
  diagnosis[RE] = exon_range[1];

  R = gene_exon_index.find(exon_range[0], exon_range[1]);

  if len(R) == 0:
    diagnosis[NE] = True;
    diagnosis[MP] = [];
    return diagnosis;
  #fi

  if len(R) > 1:
    diagnosis[IR] = True;
  #fi

  if len(R) == 1 and exon_range[0] == R[0][0] and exon_range[1] == R[0][1]:
    diagnosis[IE] = True;
  #fi

  if exon_range[0] != min([r[0] for r in R]):
    diagnosis[A5] = True;
  if exon_range[1] != max([r[1] for r in R]):
    diagnosis[A3] = True;

  diagnosis[MP] = [ r[2] for r in R ];

  return diagnosis;
#edef

###############################################################################

def diagnose_alt(alt_gene_exons, gene_exon_index, gene_exons, strand):
  diagnoses = [];
  diagnoses_w_skipped = [];
  for exon in alt_gene_exons:
    diagnoses.append(diagnose_exon(exon, gene_exon_index));
  #efor

  gene_n_exons = len(gene_exons.keys());

  exon_i = 1;
  exon_list = [ i for i in xrange(gene_n_exons) ];
  while exon_i < gene_n_exons+1:
    if len(diagnoses) > 0:
      if diagnoses[0][NE]:
        if diagnoses[0][RS] < gene_exons[exon_i][0]:
          diagnoses[0][ORS] = 0;
          diagnoses[0][ORE] = 0;
          diagnoses_w_skipped.append(diagnoses.pop(0));
          continue;
        #fi
      else:
        if exon_i in diagnoses[0][MP]:
          exon_i = exon_i + 1;
          continue;
        if exon_i > max(diagnoses[0][MP]):
          diagnoses[0][ORS] = min([ gene_exons[mpe][0] for mpe in diagnoses[0][MP] ]);
          diagnoses[0][ORE] = max([ gene_exons[mpe][1] for mpe in diagnoses[0][MP] ]);
          diagnoses_w_skipped.append(diagnoses.pop(0));
          continue;
        #fi
      #fi
    #fi
    d = [False]*len(diag_indexes);
    d[RS] = gene_exons[exon_i][0];
    d[RE] = gene_exons[exon_i][1];
    d[ORS] = gene_exons[exon_i][0];
    d[ORE] = gene_exons[exon_i][1];
    d[ES] = True;
    d[MP] = [ exon_i ];
    diagnoses_w_skipped.append(d);
    exon_i = exon_i+1;
    continue;
    #fi
  #ewhile
  while len(diagnoses) > 0:
    if not(diagnoses[0][NE]):
      diagnoses[0][ORS] = min([ gene_exons[mpe][0] for mpe in diagnoses[0][MP] ]);
      diagnoses[0][ORE] = max([ gene_exons[mpe][1] for mpe in diagnoses[0][MP] ]);
    #fi
    diagnoses_w_skipped.append(diagnoses.pop(0));
  #ewhile

    # If we are on the reverse strand
  if strand == '-':
    rev = [];
    # Fix the exon IDs (max - i + 1)

    for ex in diagnoses_w_skipped[::-1]:
      ex[MP] = [ (gene_n_exons - mpe + 1) for mpe in ex[MP][::-1] ];
      rev.append(ex);
    #efor
    diagnoses_w_skipped = rev;
  #fi

  return diagnoses_w_skipped;
#efor

###############################################################################

def print_diagnosis(D):
# 3: Intron retention
# 4: Alternative 3' SS
# 5: Alternative 5' SS
# 6: Exon skipping
# 7: New exon
# 8: Identical
  ex_d = [];
  for d in D:
    ex_d = [];
    if d[IR]:
      ex_d.append('Intron retention');
    if d[A3]:
      ex_d.append("Alternative 3' SS");
    if d[A5]:
      ex_d.append("Alternative 5' SS");
    if d[ES]:
      ex_d.append("Skipped Exon");
    if d[NE]:
      ex_d.append("New Exon");
    if d[IE]:
      ex_d.append("Identical exon");
    print d[2], ' / '.join(ex_d);
  #efor
#edef

###############################################################################

def diagnose_annotation(IDX, GN, attr_name = 'unsplitGeneID'  display=False):
  ALT_R = GN.Get(_.attr / 'geneid', *GN.Names);
  ALT_R = ALT_R.To(_.geneid, Do=_.Each(lambda x: dict([ tuple([m.strip() for m in y.split('=')]) for y in x.split(';')])[attr_name]))[_.parent != ""];
  ALT_R = ALT_R.GroupBy(_.geneid);
  ALT_R = ALT_R.Get(_.geneid[0], _.parent, _.start, _.end, _.strand).GroupBy(_.parent).To(_.strand, Do=_[0]);

  D = [];  

  for g in zip(*ALT_R()):
    n_exons, idx = IDX[g[0]];
    for alt, starts, ends, strand in zip(*g[1:]):
      alt_diag = diagnose_alt(zip(starts, ends), idx, n_exons, strand);
      for exon in alt_diag:
        D.append((g[0], alt) + tuple(exon));
      #efor

      if display:
        print g[0], alt;
        print_diagnosis(alt_diag);
      #fi
    #efor
  #efor
  return Rep(D) / ('orig_gene', 'alt_gene', 'exon_start', 'exon_end', 'orig_exon_start', 'orig_exon_end', 'orig_exon_coverage', 'retention', 'alt5', 'alt3', 'skipped', 'new_exon', 'identical');
#edef

###############################################################################

def gather_diagnosis_statistics(D):
  Da = D;
  D  = D[_.alt_gene != _.orig_gene] # Remove the original annotation, we only care about isoforms.

  S = {};

  S['n_isoforms_p_gene'] = D.Unique(_.alt_gene).GroupBy(_.orig_gene).Get(_.orig_gene, _.alt_gene.Count())
  S['n_isoforms']        = D.Unique(_.alt_gene).alt_gene.Shape();

  S['n_skipped_exon_p_exon_position'] = D[_.skipped].Flat().GroupBy(_.orig_exon_coverage).Get(_.orig_exon_coverage, _.skipped.Count());
  exon_dist_norm = Da[_.orig_gene.In(Da.Get(_.orig_gene, _.alt_gene).GroupBy(_.orig_gene)[_.alt_gene.Unique().Count() > 2][_.orig_gene == _.alt_gene].orig_gene)].Flat().GroupBy(_.orig_exon_coverage).Get(_.orig_exon_coverage, _.alt_gene.Unique().Count());
  S['n_skipped_exon_p_exon_position_norm'] = (S['n_skipped_exon_p_exon_position'] | Match | exon_dist_norm).Get(_.orig_exon_coverage, _.skipped.Cast(float) / _.alt_gene.Cast(float));
  S['length_of_skipped_exons']        = D[_.skipped].Get(_.exon_end - _.exon_start + 1);
  S['intron_retention_size']          = D[_.retention].Get(_.orig_gene, _.orig_exon_coverage.Count()).GroupBy(_.orig_exon_coverage).Get(_.orig_exon_coverage, _.orig_gene.Count());

  S['length_exon_known_isoforms']     = Da[_.orig_gene.In(Da.Get(_.orig_gene, _.alt_gene).GroupBy(_.orig_gene)[_.alt_gene.Unique().Count() > 2][_.orig_gene == _.alt_gene].orig_gene)][_.orig_gene == _.alt_gene].Get(_.orig_exon_end - _.orig_exon_start);
  S['length_exon_known_no_isoforms']  = Da[_.orig_gene == _.alt_gene].Get(_.orig_exon_end - _.orig_exon_start);
  S['length_exon_isoforms']           = D[~_.skipped].Get(_.exon_end - _.exon_start);

  S['length_alt5'] = D[_.alt5].Get(_.exon_start - _.orig_exon_start);
  S['length_alt3'] = D[_.alt3].Get(_.exon_end - _.orig_exon_end);

  S['feat_p_isoform_sum'] = gather_per_isoform_statistics(D, func=lambda d: sum(d));
  S['feat_p_isoform_any'] = gather_per_isoform_statistics(D, func=lambda d: any(d));
  S['feat_p_isoform_all'] = gather_per_isoform_statistics(D, func=lambda d: all(d));
  #S['feat_p_isoform_avg'] = gather_per_isoform_statistics(D, func=lambda d: float(sum(d))/float(len(d)));

  S['feat_p_gene_sum'] = gather_per_gene_statistics(D, func=lambda d: sum(d));
  S['feat_p_gene_any'] = gather_per_gene_statistics(D, func=lambda d: any(d));
  S['feat_p_gene_all'] = gather_per_gene_statistics(D, func=lambda d: all(d));
  #S['feat_p_gene_avg'] = gather_per_gene_statistics(D, func=lambda d: float(sum(d))/float(len(d)));

  return S;
#edef

###############################################################################

def gather_per_isoform_statistics(D, func=lambda d: any(d)):
  S = {};

  M = D.GroupBy(_.alt_gene).To(_.orig_gene, Do=_.Get(_.orig_gene[0])).Without(_.exon_start, _.exon_end, _.orig_exon_coverage);
  M = M.To(_.retention, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.alt5, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.alt3,  Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.skipped, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.new_exon, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.identical, Do=_.Array().Each(lambda x: func(x)));

  S['retention'] = M[_.retention > 0].retention.Sum()();
  S['alt5']      = M[_.alt5 > 0].alt5.Sum()();
  S['alt3']      = M[_.alt3 > 0].alt3.Sum()();
  S['skipped']   = M[_.skipped > 0].skipped.Sum()();
  S['new_exon']  = M[_.new_exon > 0].new_exon.Sum()();
  S['identical'] = M[_.identical > 0].identical.Sum()();

  return S;
#edef

###############################################################################

def gather_per_gene_statistics(D, func=lambda d: any(d)):
  S = {};

  M = D.GroupBy(_.orig_gene).Without(_.alt_gene, _.exon_start, _.exon_end, _.orig_exon_coverage);
  M = M.To(_.retention, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.alt5, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.alt3,  Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.skipped, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.new_exon, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.identical, Do=_.Array().Each(lambda x: func(x)));
  
  S['retention'] = M[_.retention > 0].retention.Sum()();
  S['alt5']      = M[_.alt5 > 0].alt5.Sum()();
  S['alt3']      = M[_.alt3 > 0].alt3.Sum()();
  S['skipped']   = M[_.skipped > 0].skipped.Sum()();
  S['new_exon']  = M[_.new_exon > 0].new_exon.Sum()();
  S['identical'] = M[_.identical > 0].identical.Sum()();

  return S;

#edef

###############################################################################

def draw_statistics(S, odir, normed_hist=True, outfmt='svg'):

  # Plot skipped exon number
    x, y = S['n_skipped_exon_p_exon_position']();
    draw_lineplot(x, y, "How often an exon at position n is skipped", "exon position", "number of times skipped", odir, outfmt=outfmt);

  # Plot probability of skipping exon n
    x, y = S['n_skipped_exon_p_exon_position_norm']();
    draw_lineplot(x, y, "Probability of skipping exon at position n", "Exon position", "Fraction skipped", odir, outfmt=outfmt);

  # Plot the length distribution of skipped exons
    x = S['length_of_skipped_exons']();
    draw_hist(x, "Length distribution of skipped exons", "Exon length", "Skipping frequency", normed_hist, odir, nbins=100, outfmt=outfmt);

  # Plot the length distribution of exons in gene with known isoforms
    x = S['length_exon_known_isoforms']();
    draw_hist(x, "Length distribution of exons in genes with known isoforms", "Exon length", "Frequency", normed_hist, odir, nbins=100, outfmt=outfmt);

  # Plot the length distribution of exons in all genes (no isoforms)
    x = S['length_exon_known_no_isoforms']();
    draw_hist(x, "Length distribution of exons in genes", "Exon length", "Frequency", normed_hist, odir, nbins=100, outfmt=outfmt);

  # Plot the length distribution of exons in isoforms (not including skipped exons)
    x = S['length_exon_isoforms']();
    draw_hist(x, "Length distribution of exons in isoforms", "Exon length", "Frequency", normed_hist, odir, nbins=100, outfmt=outfmt);

  # Plot isoform distribution (Number of isoforms per gene)
    x, y = S['n_isoforms_p_gene'].GroupBy(_.alt_gene).Get(_.alt_gene, _.orig_gene.Count())();
    draw_lineplot(x, y, "Number of isoforms per gene", "Number of isoforms", "Frequency", odir, outfmt=outfmt);

  # Distribution of retention sizes
    x, y = S['intron_retention_size']();
    draw_lineplot(x, y, 'Length (in exons) of intron retentions', 'Number of exons', 'Occurances', odir, outfmt=outfmt)

  # Distribution of the length of alternative 5' splicing sites
    x = S['length_alt5']();
    draw_hist(x, "Length of alternative 5' splicing sites", "Length (bp)", "Occurances", normed_hist, odir, 100, outfmt=outfmt);

  # Distribution of the length of alternative 3' splicing sites
    x = S['length_alt3']();
    draw_hist(x, "Length of alternative 3' splicing sites", "Length (bp)", "Occurances", normed_hist, odir, 100, outfmt=outfmt);

  # Draw table of numbers
    keys = sorted([ k for k in S.keys() if 'feat_p_' in k]);
    vals = S[keys[0]].keys();

    T = [];
    for k in keys:
      r = [ k ];
      for v in vals:
        r.append(S[k][v]);
      #efor
      T.append(tuple(r));
    #efor

    return Rep(T) / (('count',) + tuple(vals));

#edef

###############################################################################

def draw_hist(x, title="title", xlab="x", ylab="y", normed=True, odir="", nbins=None, outfmt='eps'):

  if len(x) == 0:
    return;
  #fi

  plt.cla();
  if nbins:
    plt.hist(x, nbins, normed=normed);
  else:
    plt.hist(x, normed=normed);
  #fi
  plt.xlabel(xlab);
  plt.ylabel(ylab);
  plt.title(title);

  plt.savefig('%s%s.%s' % (odir + ('/' if odir else ""), '_'.join(title.split(None)), outfmt), format=outfmt);

#edef

###############################################################################

def draw_lineplot(x, y, title="title", xlab="x", ylab="y", odir="", xlim=None, ylim=None, outfmt='eps'):

  if len(x) == 0 or len(y) == 0:
    return;
  #fi

  plt.cla();
  plt.plot(x, y, marker='x');
  plt.xlabel(xlab);
  plt.ylabel(ylab);
  plt.title(title);

  if xlim == None:
    xmin = min(x);
    xmax = max(x);
    xlim = [xmin, xmax];
  #fi

  if ylim == None:
    ymin = min(y);
    ymax = max(y);
    ylim = [ymin, ymax];
  #fi

  plt.xlim(xlim);
  plt.ylim(ylim);

  plt.savefig('%s%s.%s' % (odir + ('/' if odir else ""), '_'.join(title.split(None)), outfmt), format=outfmt);
#edef

def draw_gene_isoforms_color(retention, alt5, alt3, skipped, new, ident):
  if ident:
    return 'b';
  elif retention:
    return 'g';
  elif (alt5 or alt3):
    return 'm'
  elif new:
    return 'r';
  else:
    return 'k';
#edef

def draw_gene_isoforms(D, gene_id, odir='./', outfmt='eps'):

  import matplotlib.patches as mpatches;
  from matplotlib.collections import PatchCollection;

  ISO = D[_.orig_gene == gene_id].GroupBy(_.alt_gene).Without(_.orig_gene, _.orig_exon_start, _.orig_exon_end).Sort(_.alt_gene);

  plt.cla();

  y_loc   = 0;
  y_step  = 30;
  n_iso   = ISO.alt_gene.Shape()();
  origins = np.array([ [0, y] for y in xrange((y_step * (n_iso+1)),n_iso,-y_step) ]);
  patch_h = 10;

  xlim = [ ISO.exon_start.Min().Min()(), ISO.exon_end.Max().Max()()];
  ylim = [ y_step, (y_step * (n_iso+1)) + 2*patch_h];

  patches = [];
  
  for (origin, alt_id, starts, ends, exons, retention, alt5, alt3, skipped, new, ident) in zip(origins, *ISO()):

    plt.gca().text( min(starts), origin[1] + patch_h, alt_id, fontsize=10);
    for (exon_start, exon_end, exon_coverage, exon_retention, exon_alt5, exon_alt3, exon_skipped, exon_new, exon_ident) in zip(starts, ends, exons, retention, alt5, alt3, skipped, new, ident):
      if not(exon_skipped):
        patch = mpatches.FancyBboxPatch(origin + [ exon_start, 0], exon_end - exon_start, patch_h, boxstyle=mpatches.BoxStyle("Round", pad=0), color=draw_gene_isoforms_color(exon_retention, exon_alt5, exon_alt3, exon_skipped, exon_new, exon_ident));
        text_x, text_y = origin + [ exon_start, +patch_h/2];
        annots = zip(['Retention', "Alt 5'", "Alt 3'", "Skipped", 'New'], [exon_retention, exon_alt5, exon_alt3, exon_skipped, exon_new]);
        text  = '%s: %s' %( ','.join([str(exid) for exid in exon_coverage]), '\n'.join([ s for (s,b) in annots if b]));
        plt.gca().text(text_x, text_y, text, fontsize=10, rotation=-45);
        plt.gca().add_patch(patch);

        if all(ident):
          plt.gca().plot([exon_start, exon_start], [origin[1], 0], '--k', alpha=0.3);
          plt.gca().plot([exon_end, exon_end], [origin[1], 0], '--k', alpha=0.3);
        #fi
      #fi
    #efor
  #efor

  plt.xlim(xlim);
  plt.ylim(ylim);
  plt.title('Isoforms for gene %s' % gene_id);
  plt.xlabel('Location on chromosome');
  plt.gca().get_yaxis().set_visible(False);
  plt.gca().spines['top'].set_color('none');
  plt.gca().spines['left'].set_color('none');
  plt.gca().spines['right'].set_color('none');
  plt.tick_params(axis='x', which='both', top='off', bottom='on');
  plt.savefig('%s/isoforms_of_%s.%s' % (odir, gene_id, outfmt), format=outfmt);

  return ISO;

#edef

###############################################################################

def usage(arv0):

  print "%s: Perform isoform analysis based on two annotations"
  print "  Usage: %s <orig_gff> <new_gff> <geneid_attr> <outdir>"
  print ""
  print "    orig_gff: The original gff file, without isoform annotations"
  print "    new_gff:  The new gff file, containing isoform annotations."
  print "              Each entry should have an attribute geneID, which"
  print "              corresponds to a geneID in the original annotation."
  print "    outdir:   The directory to output figures";
  print "";
  print "  Example: %s annot.gff isoforms.gff unsplitGeneID isof_stats";

#edef

###############################################################################

def gff_add_diagnosis(D, GN):

  F = D.Get(_.orig_gene, _.alt_gene, _.exon_start, _.exon_end, _.Get(_.retention, _.alt5, _.alt3, _.skipped, _.new_exon, _.identical).Each(lambda a, b, c, d, e, f: "retention=%s;alt5=%s;alt3=%s;skipped=%s;new_exon=%s;identical=%s" % tuple([str(x) for x in [a,b,c,d,e,f]])) / 'diagnosis').Detect();
  GM = F | Match((_.alt_gene, _.exon_start, _.exon_end), (_.parent, _.start, _.end), merge_same=True, jointype='full') | GN

  FM.Get(_.seqname, _.source, _.id, _.parent, _.name, _.feature, _.start, _.end, _.score, _.strand, _.frame, (_.attr + _.diagnosis) / 'attr');

  return FM;

#edef

###############################################################################

if __name__ == 'main':

  if len(sys.argv) < 5:
    usage(sys.argv[0]);
    sys.exit(1);

  orig_gff  = sys.argv[1];
  new_gff   = sys.argv[2];
  attr_name = sys.argv[3];
  outdir    = sys.argv[4];

  GO = Read(orig_gff);
  NG = Read(new_gff);

  IDX = index_gff3_id(GO);

    # For each exon in each isoform, classify what kind of variant it is
  D = diagnose_annotation(IDX, NG, attr_name display=False);
    # Gather statistics
  S = gather_diagnosis_statistics(D);
    # Print the statistics
  draw_statistics(S, outdir, outfmt='svg');

  GD = Export(gff_add_diagnosis(D, GN), '%s/diagnosed_annotation.gff' % outdir);

  sys.exit(0);
#fi

