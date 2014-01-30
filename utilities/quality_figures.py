from pipeline_common import *
from ibidas import *
import os
import os.path
import re
import numpy
import pylab as pl
from figure_tools import *


def create_trim_figs(outdir, sample_names, filenames):
    trimlogfile = os.path.join(outdir, 'TRIMMOMATIC.std.log')
    r = getCommandOutput('cat %s | grep "Input Read"' % trimlogfile)
    r = r.split('\n')
    v = re.compile(r"Input Read Pairs\: ([\d]+) Both Surviving: ([\d]+) \([\d\.]+\%\) Forward Only Surviving: ([\d]+) \([\d\.]+\%\) Reverse Only Surviving: ([\d]+) \([\d\.]+\%\) Dropped: ([\d]+) \([\d\.]+\%\)")
    data = []
    for row in r:
        if row:
            data.append(v.match(row).groups())

    data = Rep((sample_names, data)).To(_.data, Do=_.Fields()).Detect()/('sample_names', 'input','both_surviving','forward_only','reverse_only','dropped')
    Save(data, os.path.join(outdir, 'trim_data.dat'))

    fig = dual_bargraph(data.input(), data.both_surviving(), "Number of reads", sample_names,
                  ('input read pairs','both reads surviving'), (20, 12.0))

    pl.savefig(filenames[0], dpi=200)
    #pl.show()

    ratio_surviving = (data.both_surviving.Cast(float) / data.input)()
    ratio_dropped = (data.dropped.Cast(float) / data.input)()

    fig = dual_bargraph(ratio_surviving, ratio_dropped, "Ratio of read pairs", sample_names,
                  ('surviving','dropped'), (20, 12.0))

    pl.savefig(filenames[1], dpi=200)
    #pl.show()
    
    ratio_forward = (data.forward_only.Cast(float) / data.input)()
    ratio_reverse = (data.reverse_only.Cast(float) / data.input)()

    fig = dual_bargraph(ratio_forward, ratio_reverse, "Ratio of read pairs", sample_names,
                  ('forward only surviving','reverse only surviving'), (20, 12.0))

    pl.savefig(filenames[2], dpi=200)
    #pl.show()


def create_mapping_figs(outdir, sample_names, filenames):
    umapped = re.compile(r"Uniquely mapped reads \% \|\s+([\d\.]+)\%")
    mmapped = re.compile(r"\% of reads mapped to multiple loci \|\s+([\d\.]+)\%")
    tooshort = re.compile(r"\% of reads unmapped: too short \|\s+([\d\.]+)\%")
    mismatch = re.compile(r"\% of reads unmapped\: too many mismatches \|\s+([\d\.]+)\%")
    maplength = re.compile(r"Average mapped length \|\s+([\d\.]+)")
    inputlength = re.compile(r"Average input read length \|\s+([\d\.]+)")
    data = []
    for sample_name in sample_names:
        fname = os.path.join(outdir, sample_name + ".star_align.final.log")
        r = getCommandOutput('cat %s' % fname)
        data.append(umapped.search(r).groups() + mmapped.search(r).groups() + tooshort.search(r).groups() + mismatch.search(r).groups() + maplength.search(r).groups() + inputlength.search(r).groups())

    data = Rep((sample_names, data)).To(_.data, Do=_.Fields()).Detect()/('sample_names', 'unique_match','multi_match','too_short','too_many_mismatches','avg_mapping_length', 'avg_input_length')
    print data    
    Save(data, os.path.join(outdir, 'mapping_data.dat'))

    ratio_unique = data.unique_match() / 100.0
    ratio_multi = data.multi_match() / 100.0

    fig = dual_bargraph(ratio_unique, ratio_multi, "Ratio of read pairs", sample_names,
                  ('unique match','multi match'), (20, 12.0))

    pl.savefig(filenames[0], dpi=200)
    #pl.show()

    

    ratio_short = data.too_short() / 100.0
    ratio_mismatch = data.too_many_mismatches() / 100.0

    fig = dual_bargraph(ratio_short, ratio_mismatch, "Ratio of read pairs", sample_names,
                  ('match too short','too many mismatches'), (20, 12.0))
    
    pl.savefig(filenames[1], dpi=200)
    #pl.show()


    avg_input_length = data.avg_input_length()
    avg_mapped_length = data.avg_mapping_length()
    
    fig = dual_bargraph(avg_input_length, avg_mapped_length, "Read length", sample_names,
                  ('average input length reads','average mapped length reads'), (20, 12.0))
    
    pl.savefig(filenames[2], dpi=200)
    #pl.show()

