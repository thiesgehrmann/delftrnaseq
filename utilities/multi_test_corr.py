import numpy as np;

###############################################################################

def fdr_bh(p, alpha):
  """
q = fdr_bh(p, alpha):
  Benjamini-Hochberg procedure

  p:     A list of p-values
  alpha: An alpha value

  output:
    q: A list of adjusted p-values
  """
  return fdr_gen(p, alpha, 'bh');

#edef

###############################################################################

#def fdr_bhy(p, alpha):
# """
#q = fdr_bhy(p, alpha):
#  Benjamini-Hochberg-Yekutieli procedure
#
#  p:     A list of p-values
#  alpha: An alpha value
#
#  output:
#    q: A list of adjusted p-values
#  """
#  return fdr_gen(p, alpha, 'bhy');
#
##edef

###############################################################################

def fdr_gen(p, alpha, type):

  nt     = len(p);  # Number of tests
  ps   = sorted(p);
  indx = [ i[0] for i in sorted(enumerate(p), key=lambda x:x[1]) ];

  if type == 'bhy':
    cm    = sum([ 1.0/float(i) for i in xrange(nt)] );         
    klist = [ (float(i+1)/(float(nt) * cm)) for i in xrange(nt) ];
  else:
    klist = [ (float(i+1)/float(nt)) for i in xrange(nt) ];
  #fi

    # Adjust pvalues to qvalues
  q = [ ps[i] / klist[i] for i in xrange(nt)];
    # Fix pvalues larger than 1
  q = [ qi if qi < 1.0 else 1.0 for qi in q ];

    # Monotonicity
  qm = [];
  prev_v = q[0];
  for v in q:
    qm.append(max(prev_v, v));
    prev_v = qm[-1];
  #efor

    # get back to original sorting
  qrs = [0] * nt;
  for i in xrange(nt):
    qrs[indx[i]] = qm[i];
  #efor

  return qrs;
#edef

###############################################################################
