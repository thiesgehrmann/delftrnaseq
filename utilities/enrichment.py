import numpy as np;
import scipy as sp;
from scipy import stats as ssp;
from ibidas import *;

import multi_test_corr as mtc;
reload(mtc);

###############################################################################

def fast_enrich_sample(D, M, alpha):

  #         TERM  NTERM
  #        +-----+-----+
  # DIFF  |  a  |  b  | ND
  #       +-----+-----+
  # NDIFF |  c  |  d  | NND
  #       +-----+-----+

  # D.0 = test_id
  # D.1 = significant
  # D.2 = annotation_id

  D  = D / ('test_id', 'significant', 'annotation_id');
  Df = D.FlatAll();
  Dg = Df.GroupBy(_.annotation_id);

  NG  = D.test_id.Unique().Shape()();
  ND  = D.Unique(_.test_id)[_.significant == 'yes'].test_id.Shape()();
  NND = NG - ND;
  NT  = Dg.annotation_id.Shape()();

  a = Dg.significant[_.significant == 'yes'].Shape().Get(1)();
  c = Dg.significant[_.significant == 'no'].Shape().Get(1)();

  b = [ (ND - x)  for x in a ];
  d = [ (NND - x) for x in c ];

  T = Dg.annotation_id();

  p = [ ssp.fisher_exact([ [a[i]+1,b[i]+1], [c[i]+1,d[i]+1]])[1] for i in xrange(NT)];                                                                                                                                     

    # Benjamini-Hochberg procedure
  q = mtc.fdr_bh(p, alpha);

  T = zip(T, a, b, c, d, p, q);

  R = Rep(T) / ('annotation_id', 'a', 'b', 'c', 'd', 'pvalue', 'qvalue');

  if M is not None:
    R = R | Match(0, 0, merge_same="equi") | M
  #fi

  return R.Copy();

#edef

###############################################################################

