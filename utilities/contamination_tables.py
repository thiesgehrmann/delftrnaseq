from pipeline_common import *
from ibidas import *

def make_contamination_table(outdir, sample_names, top_n = 5) :
    '''
    Produce a table with top_n contamination hits per sample that can be written into the report.
    '''
    table = []
    n_entries = len(sample_names)
    for j in range(n_entries) :
        sample = sample_names[j]
        data = Load('%s/%s.unmapped_orgs.dat' % (outdir, sample))
        entries = []
        sorted_list = [ (k[0], k[1][0]) for k in sorted(data.items(), key = lambda x: x[1][0], reverse = True) ]
        for i in xrange(min(top_n, len(sorted_list))) :
            name, count = sorted_list[i]
            name = name.strip()
            name = name.split()
            name = '\\emph{' + " ".join(name[:2]) + "} " + " ".join(name[2:])
            sample_str = sample if i == 0 else ''
            entries.append([sample_str, name, count])
        if len(entries) > 0 and j > 0 :
            table.append([''] * 3)
        for entry in entries :
            table.append(entry)
    return table
