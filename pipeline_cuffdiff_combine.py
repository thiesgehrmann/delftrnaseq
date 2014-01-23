#!/usr/bin/env python
from pipeline_common import *;
from ibidas import *
from ibidas.utils import util

C = init_conf()

def to_str(value):
    return str(value).replace('\n', ' ' )

infiles = C.__cuffdiff_output__();

files = infiles[-1]

read_group = Read(files[9]);
read_group = read_group / ('tracking_id', 'condition', 'replicate', 'raw_frags', 'internal_scaled_frags', 'external_scaled_frags', 'fpkm', 'effective_length', 'status');
read_group = read_group.To(_.tracking_id, Do=_.Cast(bytes));
read_group = read_group.To(_.condition,   Do=_.Cast(bytes));
read_group = read_group.To(_.replicate,   Do=_.Cast(bytes));
read_group = read_group.To(_.raw_frags,   Do=_.Cast(int));
read_group = read_group.To(_.internal_scaled_frags, Do=_.Cast(float));
read_group = read_group.To(_.external_scaled_frags, Do=_.Cast(float));
read_group = read_group.To(_.fpkm,        Do=_.Cast(float));
read_group = read_group.To(_.effective_length,      Do=_.Cast(bytes));


#group by tracking id, condition and replicate
read_group = read_group.GroupBy(_.tracking_id, _.condition, _.replicate, flat={0: (_.status, _.effective_length), (0,1,2):(_.raw_frags, _.internal_scaled_frags, _.external_scaled_frags, _.fpkm)}).Copy()

#first part of manual implementation of inverse HArray, needs to be added to ibidas. 
idx_all = []
name_all = []
for i in range(len(C.label_names)):
    label_name = C.label_names[i]
    names = [C.sample_names[pos] for pos, sl_idx in enumerate(C.sample_labels) if sl_idx == i]
    idx_all.extend([str(i) + '_' + str(j) for j in range(len(names))])
    name_all.extend(names)

idx_lbl = [str(i) for i in range(len(C.label_names))]
read_group = read_group.Get((_.condition + "_" + _.replicate).TakeFrom(Rep((idx_all, name_all)))/"sample_names", \
                            _.condition.TakeFrom(Rep((idx_lbl, C.label_names)))/'label_names', \
                            _.tracking_id, _.condition, _.replicate, _.raw_frags, _.external_scaled_frags)

#combine condition and replicate dimension, sort on condition first and then replicate. 
read_group = read_group.Flat().Sort((_.condition, _.replicate)).Copy()

#create separate dataset which includes data as multi-dimensional array
read_group_compact = read_group.Get(_.tracking_id/'gene_id', _.external_scaled_frags, _.sample_names, _.label_names, _.condition, _.replicate).Copy()

#second part of manual implementation of inverse HArray
names = read_group.sample_names.Each(str.lower)()
dataraw = [read_group[...,i].raw_frags/(','.join(names[i]) + '_raw') for i in range(len(names))]
datanorm = [read_group[...,i].external_scaled_frags/(names[i] + '_normalized') for i in range(len(names))]
read_group = read_group.Get(_.tracking_id/'gene_id', *(dataraw + datanorm)).Copy()

#now get gene_exp.diff data for each comparison
results = [read_group];

print files;

all_expdiff = Read(files[10]).Copy()
all_expdiff = all_expdiff / ('test_id', 'gene_id', 'gene', 'locus', 'sample_1', 'sample_2', 'status', 'value_1', 'value_2', 'log2_fold_change', 'test_stat', 'p_value', 'q_value', 'significant');
#         | -      | -      | scaffold_24:203-2284 | 0        | 1        | NOTEST   | 1.41008  | 1.09656  | -0.362794 | 0        | 1        | 1        | no
all_expdiff = all_expdiff.To(_.test_id, Do=_.Cast(bytes));
all_expdiff = all_expdiff.To(_.gene_id, Do=_.Cast(bytes));
all_expdiff = all_expdiff.To(_.gene, Do=_.Cast(bytes));
all_expdiff = all_expdiff.To(_.locus, Do=_.Cast(bytes));
all_expdiff = all_expdiff.To(_.sample_1, Do=_.Cast(int));
all_expdiff = all_expdiff.To(_.sample_2, Do=_.Cast(int));
all_expdiff = all_expdiff.To(_.status, Do=_.Cast(bytes));
all_expdiff = all_expdiff.To(_.value_1, Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.value_2, Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.log2_fold_change, Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.test_stat, Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.p_value, Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.q_value, Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.significant, Do=_.Cast(bytes));

for i in xrange(len(C.cuffdiff_cmp)):
  a, b = C.cuffdiff_cmp[i];
  #files = infiles[i]
  cname = (C.label_names[a] + "_" + C.label_names[b]).lower()
  expdiff = all_expdiff[(_.sample_1==a) & (_.sample_2==b)].Copy()

  #expdiff = Read(files[6])
  expdiffd = expdiff.Without(_.gene_id, _.test_id, _.locus, _.gene, _.sample_1, _.sample_2)
  expdiffd = expdiffd / tuple([cname + '_' + name for name in expdiffd.Names])
  res = expdiff.Get(_.gene_id, expdiffd)

  results.append(res.Copy());

#efor

#combine the data
for i in range(1, len(results)):
    results[0] = (results[0] |Match(_.gene_id, _.gene_id)| results[i]).Copy()
#efor

#if annots file is available, add that one too (should have gene_id slice). 
if C.annots_file:
    annots = Load(C.annots_file)
    res = annots |Match| results[0]
else:
    res = results[0]
res = res % 'genes'

res2 = res |Match| read_group_compact
res2 = res2 % 'genes'

Save(res2, C.__cuffdiff_combine_output__()[0])
Save(res, C.__cuffdiff_combine_output__()[1])

sys.exit(0);
