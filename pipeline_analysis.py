#!/usr/bin/env python

import sys, os

from pipeline_common import *;
from ibidas import *
from ibidas.utils import util

sys.path.append('%s/utilities/' % os.path.dirname(os.path.realpath(__file__)));
import analysis_figures as af
import latex

C = init_conf()

###############################################################################


infiles = C.__cuffdiff_combine_output__();
data    = Load(infiles[0])

outfiles = C.__analysis_output__()

for (i, filter) in enumerate(C.analysis_filter):

  print "Performing analysis for filter %s." % (str(filter));

  filt_outfiles = outfiles[i];

  if filter == None:
    data_f = data;
  else:
    filt = Read(filter);
    data_f = data[_.test_id.In(filt.Get(0))];
  #fi

  l = latex.LatexFile(filt_outfiles[-1][0])
  l.write_title('Results "%s"' % C.title, C.author)
  l.add_text("Filter for this analysis: %s.\n" % str(filter));

  l.start_section("Significant genes")
  difffile = filt_outfiles[2][0]
  af.create_diffgenes_stats(data_f, difffile)
  l.include_figure(difffile, 'diffgenes', 'Number of significant genes for each comparison. Upregulated/Downregulated means that the gene has a respectively higher/lower expression in the second mentioned condtion, compared to the first mentioned condition.', width=1.5)

  venn_subsets_file = filt_outfiles[-2][0];
  venn_files        = filt_outfiles[0];
  subsets = [];
  for pos, (venn, filenames) in enumerate(zip(C.analysis_venn, venn_files)):
    ss = af.create_venn(data_f, [C.cuffdiff_cmp[v] for v in venn], C.label_names, filenames, C.analysis_venn_updown_split)
    subsets.extend(ss);
    for filename in filenames:
      l.include_figure(filename, 'venn%d' % pos, "Overlap significant genes between condition comparisons.")
    #efor
  #efor
  Save(Rep(subsets) / ('group', 'subset', 'genes'), venn_subsets_file);
  l.clear_page()
  l.start_section("Patterns across all genes")

  pcafile = filt_outfiles[1][0]
  af.create_pca(data_f, pcafile)
  l.include_figure(pcafile, 'pca', 'Distances between conditions after PCA dimension reduction (applied to normalized frag count data).')

  pcafile2 = filt_outfiles[1][1]
  af.create_pca(data_f, pcafile2, logfilter=True)
  l.include_figure(pcafile2, 'pcalog', 'Distances between conditions after PCA dimension reduction (applied to log transformed, filtered, normalized frag count data). Frag counts were transformed to log space using $log2(fragcount + 16)$. Genes were filtered on having a $std > 0.5$. ')

  clustersfile = filt_outfiles[3][0]
  af.create_cluster(data_f, clustersfile)
  l.include_figure(clustersfile, 'clustering', 'k-Means clustering (n=12 clusters) applied to all genes. Normalized frag counts per gene per replicate were transformed to log space using log2(fragcount + 16). Genes were filtered on having a $std > 0.5$. Remaining genes were standardized (mean = 0, std = 1), and subsequently clustered.',width=2.0)

  l.end_document()
#efor

sys.exit(0);
