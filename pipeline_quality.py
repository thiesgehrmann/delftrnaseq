#!/usr/bin/env python

import sys;

from pipeline_common import *;
from ibidas import *
from ibidas.utils import util

sys.path.append('%s/utilities/' % os.path.dirname(os.path.realpath(__file__)));
import quality_figures as qf
import latex


C = init_conf()

outfiles = C.__quality_output__()

l = latex.LatexFile(outfiles[-1][0])
l.write_title('Quality Report "%s"' % C.title, C.author)

l.start_section('Read Filtering')

qf.create_trim_figs(C.outdir, C.sample_names, outfiles[0])
l.include_figure(outfiles[0][0], 'nreads', 'Number of read pairs before and after filtering using Trimmomatic. Adapters and low-quality regions in the reads are removed. Reads that become too short due to the filtering are dropped.')
l.include_figure(outfiles[0][1], 'ratio', 'Ratio of read pairs were both reads survive filtering through Trimmomatic, and ratio of read pairs for which both reads are dropped.')
l.include_figure(outfiles[0][2], 'single', 'Ratio of read pairs for which only the forward or only the reverse read survives filtering through Trimmomatic.')

l.clear_page()

l.start_section("Read Mapping")
qf.create_mapping_figs(C.outdir, C.sample_names, outfiles[1])
l.include_figure(outfiles[1][0], 'mapped', 'Ratio of reads that is mapped uniquely or non-uniquely to the target genome.')
l.include_figure(outfiles[1][1], 'unmapped', 'Ratio of reads that is dropped due to the found alignment being too short, or due there being too many mismatches in the alignment.')
l.include_figure(outfiles[1][2], 'length', 'The average input read length compared to the average length of these reads that is aligned to the genome.')

l.clear_page()

if C.check_contamination:
  l.start_section("Read contamination")
#fi

l.end_document()

sys.exit(0);
