from pipeline_common import *
from Bio import SeqIO
import latex
import os
from ibidas import *
import re

def shorten_timmomatic_opts(opts) :
    args = opts.split(':')
    args[1] = os.path.basename(args[1])
    return ':'.join(args)

def make_adapter_table(trimmomatic_opts) :
    '''
    Produce a table with trimmed adapters.
    '''
    args = trimmomatic_opts.split(':')
    filename = args[1]
    table = []
    fin = open(filename, 'rU')
    for record in SeqIO.parse(fin, 'fasta') :
        name, seq = record.name, record.seq.tostring()
        table.append([name, '\\verb=%s=' % seq])
    fin.close()
    return table

def make_read_count_table(outdir, sample_names) :
    mapped_regex = re.compile(r"Uniquely mapped reads number \|\s+([\d\.]+)")
    multimapped_regex = re.compile(r"Number of reads mapped to multiple loci \|\s+([\d\.]+)")
    data = []
    for sample_name in sample_names:
        fname = os.path.join(outdir, sample_name + ".star_align.final.log")
        r = getCommandOutput('cat %s' % fname)
        fname = os.path.join(outdir, sample_name + '.star_align_sort.count')
        cds_count = getCommandOutput("cat  %s |  grep -v '__no_feature' | grep -v '__not_aligned' | grep -v '__alignment_not_unique' | grep -v '__too_low_aQual' | awk 'BEGIN {s=0} {s=s+$2} END {print s}'" % fname)
        data.append(mapped_regex.search(r).groups() + multimapped_regex.search(r).groups() + (int(cds_count), ))
    data = Rep((sample_names, data)).To(_.data, Do=_.Fields()).Detect()/('sample_names', 'unique_match','multi_match', 'cds_align')
    
    datafile = os.path.join(outdir, 'trim_data.dat')
    data_trim = Load(datafile)
    data = data_trim.Match(data, _.sample_names, _.sample_names)
    table = []
    print data
    for sample_name in sample_names :
        rows = data.Filter(_.sample_names == sample_name).Tuple()()
        row = rows[0]
        name, input_pairs, survived_pairs, aligned, multi_aligned, cds_aligned = row[0], row[1], row[2], row[6], row[7], row[8]
        aligned += multi_aligned
        table_row = [name, '$%s$' % '{:,.0f}'.format(input_pairs).replace(',', ',\\!'), '$%s$ ($%.2f$\\%%)' % ('{:,.0f}'.format(survived_pairs).replace(',', ',\\!'), survived_pairs / float(input_pairs) * 100.0), '$%s$ ($%.2f$\\%%)' % ('{:,.0f}'.format(aligned).replace(',', ',\\!'), aligned / float(input_pairs) * 100.0), '$%s$ ($%.2f$\\%%)' % ('{:,.0f}'.format(cds_aligned).replace(',', ',\\!'), cds_aligned / float(input_pairs) * 100.0)]
        table.append(table_row)
    return table
