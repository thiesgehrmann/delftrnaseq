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

qf.create_trim_figs(C.outdir, C.sample_names, outfiles[0], C.PE)
qf.create_mapping_figs(C.outdir, C.sample_names, outfiles[1])

l = latex.LatexFile(outfiles[-1][0])
l.write_title("Quality Report ``%s''" % C.title, C.author)

l.start_section('Read Filtering')
l.add_text('Trimmomatic options:\\\\ \n\\verb=%s=' % misc.shorten_timmomatic_opts(C.trimmomatic_trim), texcape = False)
table = misc.make_adapter_table(C.trimmomatic_trim)
l.write_table(['Primer', 'Sequence'], table, "Filtered adapter sequences.", alignment = "ll")

table = misc.make_read_count_table(C.outdir, C.sample_names)
l.write_table(['Read set', 'Read pairs', 'Survived pairs', 'Aligned pairs', 'CDS-aligned pairs'], table, "Read set statistics.", alignment = "lrrrr")

for sample in C.sample_names :
    figure_filename = sample + '.star_align_sort.read_stats'
    l.include_figure(figure_filename, 'fig:read-stats-%s' % sample, 'Read (fragment) length distribution for %s sample.' % sample, additional_options = 'type=pdf,ext=.pdf,read=.pdf')

l.include_figure(outfiles[0][0], 'fig:nreads', 'Number of read pairs before and after filtering using Trimmomatic. Adapters and low-quality regions in the reads are removed. Reads that become too short after trimming are dropped.')
l.include_figure(outfiles[0][1], 'fig:ratio', 'Ratio of read pairs were both reads survive filtering using  Trimmomatic, and ratio of read pairs for which both reads are dropped.')
if C.PE :
  l.include_figure(outfiles[0][2], 'fig:single', 'Ratio of read pairs for which only the forward or only the reverse read survives filtering through Trimmomatic.')
#fi

l.clear_page()

l.start_section("Read Mapping")
l.include_figure(outfiles[1][0], 'fig:mapped', 'Ratio of reads that are mapped uniquely or non-uniquely to the target genome.')
l.include_figure(outfiles[1][1], 'fig:unmapped', 'Ratio of reads that is dropped due to the found alignment being too short, or due to having too many mismatches in the alignment.')
l.include_figure(outfiles[1][2], 'fig:length', 'The average input read length compared to the average length of these reads that is aligned to the genome.')

l.clear_page()

if C.check_contamination :
    l.start_section("Contamination")
    l.add_text('Reads that were not mapped by STAR were assembled using the Trinity assembler into transcripts. Predicted ORF on these transcripts were then BLASTed against known sequences to and the source organism of the hits was determined.')
    table = ct.make_contamination_table(C.outdir, C.sample_names, top_n = 5)
    l.write_table(['Sample', 'Organism', 'Hits'], table, "Top 5 contaminants for each sample.", alignment = "lll")
    l.clear_page()

l.start_section("FastQC: mapped reads")
for sample in C.sample_names :
    l.start_section('Sample %s' % sample, level = 1)
    dirname = C.outdir + '/' + sample + '.star_align_fastqc'
    l.include_figure(dirname + '/Images/per_base_quality.png', 'fig:' + dirname + '-quality', "Per-base read quality for %s" % sample, width = 0.6)
    l.include_figure(dirname + '/Images/kmer_profiles.png', 'fig:' + dirname + '-kmers', "$k$-mer over-representation profiles for %s" % sample, width = 0.6)
    l.clear_page()

l.start_section("FastQC: unmapped reads")
for sample in C.sample_names :
    l.start_section('Sample %s' % sample, level = 1)
    if C.PE :
        dirname = C.outdir + '/' + sample + '.star_align_unmapped_R1_fastqc'
        l.include_figure(dirname + '/Images/per_base_quality.png', 'fig:' + dirname + '-quality', "Per-base read quality for %s (read 1)" % sample, width = 0.6)
        l.include_figure(dirname + '/Images/kmer_profiles.png', 'fig:' + dirname + '-kmers', "$k$-mer over-representation profiles for %s (read 1)" % sample, width = 0.6)
        dirname = C.outdir + '/' + sample + '.star_align_unmapped_R2_fastqc'
        l.include_figure(dirname + '/Images/per_base_quality.png', 'fig:' + dirname + '-quality', "Per-base read quality for %s (read 2)" % sample, width = 0.6)
        l.include_figure(dirname + '/Images/kmer_profiles.png', 'fig:' + dirname + '-kmers', "$k$-mer over-representation profiles for %s (read 2)" % sample, width = 0.6)
    else :
        dirname = C.outdir + '/' + sample + '.star_align_unmapped_fastqc'
        l.include_figure(dirname + '/Images/per_base_quality.png', 'fig:' + dirname + '-quality', "Per-base read quality for %s" % sample, width = 0.6)
        l.include_figure(dirname + '/Images/kmer_profiles.png', 'fig:' + dirname + '-kmers', "$k$-mer over-representation profiles for %s" % sample, width = 0.6)
    l.clear_page()


l.end_document()


print "QUALITY REPORT is at %s. Run pdflatex twice." % (outfiles[-1][0]);

sys.exit(0);
