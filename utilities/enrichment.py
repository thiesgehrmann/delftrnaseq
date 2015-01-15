import numpy as np;
import scipy as sp;
from scipy import stats as ssp;
from ibidas import *;
from ibidas.utils.util import debug_here
import numpy;
import sys;

import multi_test_corr as mtc;
reload(mtc);

###############################################################################

def fast_enrich(D, M, alpha):

    #      any   bool      any(array)
  D = D / ('id', 'class1', 'class2');
  D = D.To(_.class1, Do=_.Cast(bool));
  D = D.ReplaceMissing();

  Dflat        = D.FlatAll();
  Dgroupclass2 = Dflat.GroupBy(_.class2);

  Nids    = D.id.Unique().Shape()();
  Nclass2 = Dgroupclass2.class2.Shape()();

  Ntrue  = D[_.class1].id.Shape()();
  Nfalse = Nids - Ntrue;

  a = Dgroupclass2[_.class1].id.Shape().Get(1)();
  c = Dgroupclass2[~_.class1].id.Shape().Get(1)();
  b = [ (Ntrue - x)  for x in a ];
  d = [ (Nfalse - x) for x in c ];

  sys.stdout.flush()

  print "Ntrue: ", Ntrue, " | Nfalse: ", Nfalse;
  print "Nids: ", Nids;

  if any(numpy.concatenate((a,b,c,d)) < 0):
    print "ERROR, SOMETHING WEIRD I CAN't EXPLAIN YET!";
    p = [ 1.0 for i in xrange(Nclass2) ];
  else:
    p = [ ssp.fisher_exact([ [a[i]+1,b[i]+1], [c[i]+1,d[i]+1]])[1] for i in xrange(Nclass2)];
  #fi
    
    # Benjamini-Hochberg procedure
  q = mtc.fdr_bh(p, alpha);

  T = zip(Dgroupclass2.class2(), a, b, c, d, p, q);

  R = Rep(T) / ('annotation_id', 'a', 'b', 'c', 'd', 'pvalue', 'qvalue');

  if M is not None:
    R = R | Match(0, 0, merge_same="equi") | M;
  #fi

  return R.Copy();
#edef

###############################################################################

def fast_enrich_sample(D, M, alpha, all_or_up_or_down='all'):
  #         TERM  NTERM
  #        +-----+-----+
  # DIFF  |  a  |  b  | ND
  #       +-----+-----+
  # NDIFF |  c  |  d  | NND
  #       +-----+-----+

  # D.0 = test_id
  # D.1 = significant
  # D.2 = annotation_id
  output_slice_names = ('annotation_id', 'a', 'b', 'c', 'd', 'pvalue', 'qvalue');

  D = D / ('test_id', 'significant', 'logfold', 'annotation_id');
  D = D.To(_.significant, Do=_.Cast('bytes'));
  D = D.To(_.logfold, Do=_.Cast('real64'));
  D = D.To(_.significant, _.logfold, Do=_.ReplaceMissing());

  if D.annotation_id.Shape().Get(1).Sum()() == 0:
    S = Rep(tuple([0 for x in output_slice_names])) / output_slice_names;
    return S[_.test_id > 0];
  #fi

    # The data grouped by annotation ids
  Df = D.FlatAll();
  Dg = Df.GroupBy(_.annotation_id);

    # The total number of genes and terms 
  NG  = D.test_id.Unique().Shape()();
  NT  = Dg.annotation_id.Shape()();

  # ND: The number of significant test_ids
  # a:  The set of significant genes
  # c:  The set of non-significant genes
  if all_or_up_or_down == 'all':
      # ND: The number of significant test_ids (differentially expressed) ND = Number Diff
      # a:  The set of significant genes
      # c:  The set of non-significant genes
    ND = D[_.significant == 'yes'].test_id.Shape()();
    a  = Dg[  _.significant == 'yes' ];
    c  = Dg[~(_.significant == 'yes')];
  elif all_or_up_or_down == 'up':
      # ND: The number of significant test_ids which are upregulated
      # a:  The set of significant genes which are upregulated
      # c:  The set of genes which are non-significant or are not up-regulated
    ND = D[  ( _.significant == 'yes' ) & ( _.logfold > 0 ) ].test_id.Shape()();
    a  = Dg[ ( _.significant == 'yes' ) & ( _.logfold > 0 ) ];
    c  = Dg[~(( _.significant == 'yes' ) & ( _.logfold > 0 )) ];
  elif all_or_up_or_down == 'down':
      # ND: The number of significant test_ids which are downregulated
      # a:  The set of significant genes which are downregulated
      # c:  The set of genes which are non-significant or are not downregulated
    ND = D[  ( _.significant == 'yes' ) & ( _.logfold < 0 ) ].test_id.Shape()();
    a  = Dg[ ( _.significant == 'yes' ) & ( _.logfold < 0 ) ];
    c  = Dg[~(( _.significant == 'yes' ) & ( _.logfold < 0 )) ];
  #fi

    # The number of genes which are not significant and not up or downregulated
  NND = NG - ND;

  a = a.test_id.Shape().Get(1)();
  c = c.test_id.Shape().Get(1)();

  b = [ (ND - x)  for x in a ];
  d = [ (NND - x) for x in c ];

  T = Dg.annotation_id();

  p = [ ssp.fisher_exact([ [a[i]+1,b[i]+1], [c[i]+1,d[i]+1]])[1] for i in xrange(NT)];                                                                                                                                     

    # Benjamini-Hochberg procedure
  q = mtc.fdr_bh(p, alpha);

  T = zip(T, a, b, c, d, p, q);

  R = Rep(T) / ('annotation_id', 'a', 'b', 'c', 'd', 'pvalue', 'qvalue');

  if M is not None:
    R = R | Match(0, 0, merge_same="equi") | M;
  #fi

  return R.Copy();

#edef

###############################################################################

