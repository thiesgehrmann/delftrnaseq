#!/usr/bin/env python

import os;

from pipeline_common import *;
from ibidas import *
from ibidas.utils import util

os.path.append('utilities');
import analysis_figures as af
import latex


C = init_conf()

###############################################################################

infiles = C.__cuffdiff_combine_output__();
data = Load(infiles[0])

outfiles = C.__analysis_output__()

l = latex.LatexFile(outfiles[-1][0])
l.write_title('Results "%s"' % C.title, C.author)

l.start_section("Significant genes")
difffile = outfiles[2][0]
af.create_diffgenes_stats(data, difffile)
l.include_figure(difffile, 'diffgenes', 'Number of significant genes for each comparison. Upregulated/Downregulated means that the gene has a respectively higher/lower expression in the second mentioned condtion, compared to the first mentioned condition.', width=1.5)


venn_files = outfiles[0]
for pos, (venn, filename) in enumerate(zip(C.analysis_venn, venn_files)):
    af.create_venn(data, [C.cuffdiff_cmp[v] for v in venn], filename)
    l.include_figure(filename, 'venn%d' % pos, "Overlap significant genes between condition comparisons.")

l.clear_page()
l.start_section("Patterns across all genes")

pcafile = outfiles[1][0]
af.create_pca(data, pcafile)
l.include_figure(pcafile, 'pca', 'Distances between conditions after PCA dimension reduction (applied to normalized frag count data).')

pcafile2 = outfiles[1][1]
af.create_pca(data, pcafile2, logfilter=True)
l.include_figure(pcafile2, 'pcalog', 'Distances between conditions after PCA dimension reduction (applied to log transformed, filtered, normalized frag count data). Frag counts were transformed to log space using $log2(fragcount + 16)$. Genes were filtered on having a $std > 0.5$. ')

clustersfile = outfiles[3][0]
af.create_cluster(data, clustersfile)
l.include_figure(clustersfile, 'clustering', 'k-Means clustering (n=12 clusters) applied to all genes. Normalized frag counts per gene per replicate were transformed to log space using log2(fragcount + 16). Genes were filtered on having a $std > 0.5$. Remaining genes were standardized (mean = 0, std = 1), and subsequently clustered.',width=2.0)

l.end_document()

sys.exit(0);
