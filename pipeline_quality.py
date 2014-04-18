#!/usr/bin/env python

import sys;

from pipeline_common import *;
from ibidas import *
from ibidas.utils import util

sys.path.append('%s/utilities/' % os.path.dirname(os.path.realpath(__file__)));
import quality_figures as qf
import contamination_tables as ct
import misc
import latex


C = init_conf()

outfiles = C.__quality_output__()

l = latex.LatexFile(outfiles[-1][0])
l.write_title("Quality Report ``%s''" % C.title, C.author)

l.start_section('Read Filtering')
l.add_text('Trimmomatic options:\\\\ \n\\verb=%s=' % misc.shorten_timmomatic_opts(C.trimmomatic_trim))
table = misc.make_adapter_table(C.trimmomatic_trim)
l.write_table(['Primer', 'Sequence'], table, "Filtered adapter sequences.", alignment = "ll")

table = misc.make_read_count_table(C.outdir, C.sample_names)
l.write_table(['Read set', 'Read pairs', 'Survived pairs', 'Aligned pairs', 'CDS-aligned reads'], table, "Read set statistics.", alignment = "lrrrr")

qf.create_trim_figs(C.outdir, C.sample_names, outfiles[0], C.PE)
l.include_figure(outfiles[0][0], 'nreads', 'Number of read pairs before and after filtering using Trimmomatic. Adapters and low-quality regions in the reads are removed. Reads that become too short after trimming are dropped.')
l.include_figure(outfiles[0][1], 'ratio', 'Ratio of read pairs were both reads survive filtering using  Trimmomatic, and ratio of read pairs for which both reads are dropped.')
if C.PE :
  l.include_figure(outfiles[0][2], 'single', 'Ratio of read pairs for which only the forward or only the reverse read survives filtering through Trimmomatic.')
#fi

l.clear_page()

l.start_section("Read Mapping")
qf.create_mapping_figs(C.outdir, C.sample_names, outfiles[1])
l.include_figure(outfiles[1][0], 'mapped', 'Ratio of reads that are mapped uniquely or non-uniquely to the target genome.')
l.include_figure(outfiles[1][1], 'unmapped', 'Ratio of reads that is dropped due to the found alignment being too short, or due to having too many mismatches in the alignment.')
l.include_figure(outfiles[1][2], 'length', 'The average input read length compared to the average length of these reads that is aligned to the genome.')

l.clear_page()

if C.check_contamination :
    l.start_section("Contamination")
    l.add_text('Reads that were not mapped by STAR were assembled using the Trinity assembler into transcripts. Predicted ORF on these transcripts were then BLASTed against known sequences to and the source organism of the hits was determined.')
    table = ct.make_contamination_table(C.outdir, C.sample_names, top_n = 5)
    l.write_table(['Sample', 'Organism', 'Hits'], table, "Top 5 contaminants for each sample.", alignment = "lll")
    l.clear_page()

l.end_document()


print "QUALITY REPORT is at %s. Run pdflatex twice." % (outfiles[-1][0]);

sys.exit(0);
