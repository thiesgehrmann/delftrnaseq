#!/home/nfs/thiesgehrmann/src/ibidas/run

import sys;
from ibidas import *;

#NF, NG = split_genome(sys.argv[1], sys.argv[2]);
#Export(NF, tdir + '/genome.fasta');
#Export(NG, tdir + '/annotation.gff');
#US = unsplit(sys.argv[1], tdir + '/annotation.gff')

def split_genome(gff_file, fasta_file):

  #gff_file   = "/data/tmp/thiesgehrmann/rnaseq_schco3/SCHCO_RNASEQ.cleaned.gff";
  #fasta_file = "/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/genome.fa";
  
  G = Read(gff_file);
  F = Read(fasta_file);

    # Keep the original locations of the genes
  G = G.To(_.attr, Do=_.Get(_.attr, G.seqname, G.start, G.end).Each(lambda w, x, y, z: ';'.join(['original_seqname=%s;original_start=%d;original_end=%d' % (x, y, z), w])).Cast(bytes) / 'attr');

  GE = G[_.feature!='mRNA'].GroupBy(_.parent).Get(_.seqname[0], _.source[0], _.id[0], _.parent, _.name[0], _.feature, _.start, _.end, _.score[0], _.strand[0], _.frame, _.attr)
  
  GEF = GE | Match(_.seqname, _.f0, merge_same=True) | F;
  GEFS = GEF.To(_.seq, Do=_.Get(GEF.start.Min(), GEF.end.Max(), _.seq).Each(lambda x,y,z: z[x-1:y]) / ('seq'));
  
  GEFSR = GEFS.To(_.start, Do=_.Get(_.start - Min(_.start) + 1) / ('start')).To(_.end, Do=_.Get(_.end - GEF.start.Min() + 1) / ('end'))
  
  NF = GEFSR.Get(_.parent, _.seq).To(_.seq, Do=_.Cast('DNA'));
  NG = GEFSR.Get(_.parent / ('seqname'), _.source, _.id, _.parent, _.name, _.feature, _.start, _.end, _.score, _.strand, _.frame, _.attr).Flat();

  return NG, NF;
#edef

def split_genome_smart(gff_file, fasta_file, readlen, IS_mean, IS_stdev):

  G = Read(gff_file);
  F = Read(fasta_file);
    
    # Keep the original locations of the genes
  G = G.To(_.attr, Do=_.Get(_.attr, G.seqname, G.start, G.end).Each(lambda w, x, y, z: ';'.join(['original_seqname=%s;original_start=%d;original_end=%d' % (x, y, z), w])).Cast(bytes) / 'attr');

  neighbors    = G[_.parent != ""].GroupBy(_.parent).Get(_.seqname[0], _.parent, _.start.Min(), _.end.Max()).Sort(_.seqname, _.start).GroupBy(_.seqname);
  N            = [];
  padding_size = int(readlen*2 + IS_mean + IS_stdev -1);

  for chromosome, parents, starts, ends in zip(*neighbors()):
    K = [ ('', 0, 0) ] + zip(parents, starts, ends) ;
    K = K + [ ('', K[-1][2] + padding_size, 0) ];

    i = 0
    for k in xrange(1, len(K) - 1):
      left_dist  = K[k][1] - K[k-1][2];
      right_dist = K[k+1][1] - K[k][2];
      n = (K[k][0], min(max(int(left_dist / 2), 0), padding_size), min(max(int(right_dist / 2), 0), padding_size));
      N.append(n);
    #efor
  #efor
  N = Rep(N) / ('parent', 'left_padding', 'right_padding');

  GE = G[_.feature!='mRNA'].GroupBy(_.parent).Get(_.seqname[0], _.source[0], _.id[0], _.parent, _.name[0], _.feature, _.start, _.end, _.score[0], _.strand[0], _.frame, _.attr)
  
  GE = GE | Match(_.parent, _.parent) | N;

  GEF = GE | Match(_.seqname, _.f0, merge_same=True) | F;
  GEFS = GEF.To(_.seq, Do=_.Get(GEF.start.Min() - GEF.left_padding, GEF.end.Max() + GEF.right_padding, _.seq).Each(lambda x,y,z: z[x-1:y]) / ('seq'));
  
  GEFSR = GEFS.To(_.start, Do=_.Get(_.start - Min(_.start) + 1) / ('start')).To(_.end, Do=_.Get(_.end - GEF.start.Min() + 1) / ('end'))
  
  NF = GEFSR.Get(_.parent, _.left_padding, _.right_padding, _.seq).To(_.seq, Do=_.Cast('DNA'));
  NG = GEFSR.Get(_.parent / ('seqname'), _.source, _.id, _.parent, _.name, _.feature, (_.start + _.left_padding) / 'start', (_.end + _.left_padding) / 'end', _.score, _.strand, _.frame, _.attr).Flat();

  return NG.Copy(), NF.Copy();
#edef


def unsplit_smart(gff_file, fasta_file, created_gff):

  G  = Read(gff_file);
  NF = Read(fasta_file, sep=['|']) / ('parent', 'left_padding', 'right_padding', 'seq');
  NF = NF.To(_.left_padding, Do=_.Cast(int)).To(_.right_padding, Do=_.Cast(int));
  NG = Read(created_gff);

  CNG = NG | Match(_.seqname, _.parent2) | NF.Get(_.parent / 'parent2', _.left_padding);
  CNG = CNG.To(_.start, Do=_.Get(CNG.start, CNG.left_padding).Each(lambda x, y: x - y) / 'start');
  CNG = CNG.To(_.end,   Do=_.Get(CNG.end, CNG.left_padding).Each(lambda x, y: x - y) / 'end');
  CNG = CNG.Without(_.left_padding, _.parent2);

  NI = G[_.id != ''].Get(_.seqname, _.id, _.start) / ('oseq', 'oid', 'ostart');

  X = CNG | Match(_.seqname, _.oid) | NI;
  X = X.Get(_.oseq, _.source, _.id, _.parent, _.name, _.feature, (_.start + _.ostart - 1), (_.end + _.ostart - 1), _.score, _.strand, _.frame, _.attr);

  return X / tuple(G.Names);
#edef

def unsplit(gff_file, created_gff):

  G = Read(gff_file);
  N = Read(created_gff);

  NI = G[_.id != ''].Get(_.seqname, _.id, _.start) / ('oseq', 'oid', 'ostart');

  X = N | Match(_.seqname, _.oid) | NI;
  X = X.Get(_.oseq, _.source, _.id, _.parent, _.name, _.feature, (_.start + _.ostart - 1), (_.end + _.ostart - 1), _.score, _.strand, _.frame, _.attr);

  return X / tuple(G.Names);
#edef

def set_x2y(gff_file, x, xtype, y, ytype):
  G = Read(gff_file);

  if xtype == 'attribute' and ytype == 'name' and y in G.Names:
    print 'lol'
    G = G.To(_.attr, Do=_.Get(G.attr, G.Get(y)).Each(lambda a, t: ';'.join(['='.join(m) for m in dict([ tuple(ai.split('=')) for ai in a.split(';')] + [(x, str(t))]).items()])) / 'attr');
  else:
    return None;
  #fi
  #if xtype == 'name' and ytpe = 'name' and x in G.Names and y in G.Names:

  return G;
#edef


def usage():
  print "Transform the genome and manipulate attributes";
  print "%s split   <gff_file> <genome_fasta_file> <IS_mean> <IS_stdev> <split_output_gff> <split_output_fasta>" % sys.argv[0];
  print "%s unsplit <gff_file> <split_fasta> <split_gff_file> <unsplit_output_gff>" % sys.argv[0];
  print "%s setx2y  <split_gff_file> <x> <xtype name|attribute> <y> <ytype name|attribute> <output_gff>" % sys.argv[0];

  print '\n'.join(sys.argv);
#edef

###############################################################################

if __name__ == '__main__':
  if sys.argv[1] == 'split' and (len(sys.argv) == 9):
    print sys.argv[2], sys.argv[3];
    G,F = split_genome_smart(sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]));
    Export(G, sys.argv[7]);
    Export(F, sys.argv[8]);
  elif sys.argv[1] == 'unsplit' and (len(sys.argv) == 6):
    G = unsplit_smart(sys.argv[2], sys.argv[3], sys.argv[4]);
    Export(G, sys.argv[5]);
  elif sys.argv[1] == 'setx2y' and (len(sys.argv) == 8):
    G = set_x2y(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]);
    Export(G, sys.argv[7]);
  else:
    usage();
    sys.exit(1);
  #fi
  sys.exit(0);
#fi
  
