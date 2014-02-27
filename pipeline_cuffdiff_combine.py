#!/usr/bin/env python
from pipeline_common import *;
from ibidas import *
from ibidas.utils import util

C = init_conf()

###############################################################################

files = C.__cuffdiff_output__();

read_group_index = C.__cuffdiff_test_group_index__[C.cuffdiff_test_type];
diff_data_index  = C.__cuffdiff_test_diff_index__[C.cuffdiff_test_type];
gslice           = 'gene_id' if C.cuffdiff_test_type == 'gene' else 'test_id' # grouping_slice 


read_group = Read(files[read_group_index]);
read_group = read_group / ('tracking_id', 'condition', 'replicate', 'raw_frags', 'internal_scaled_frags', 'external_scaled_frags', 'fpkm', 'effective_length', 'status');
read_group = read_group.To(_.tracking_id,           Do=_.Cast(bytes));
read_group = read_group.To(_.condition,             Do=_.Cast(bytes));
read_group = read_group.To(_.replicate,             Do=_.Cast(bytes));
read_group = read_group.To(_.raw_frags,             Do=_.Cast(float));
read_group = read_group.To(_.internal_scaled_frags, Do=_.Cast(float));
read_group = read_group.To(_.external_scaled_frags, Do=_.Cast(float));
read_group = read_group.To(_.fpkm,                  Do=_.Cast(float));
read_group = read_group.To(_.effective_length,      Do=_.Cast(bytes));
read_group = read_group.Copy();

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
#efor

idx_lbl = [str(i) for i in range(len(C.label_names))]
read_group = read_group.Get((_.condition + "_" + _.replicate).TakeFrom(Rep((idx_all, name_all)))/"sample_names", \
                             _.condition.TakeFrom(Rep((idx_lbl, C.label_names)))/'label_names', \
                             _.tracking_id, _.condition, _.replicate, _.raw_frags, _.external_scaled_frags)

#combine condition and replicate dimension, sort on condition first and then replicate. 
read_group = read_group.Flat().Sort((_.condition, _.replicate)).Copy()

#create separate dataset which includes data as multi-dimensional array
read_group_compact = read_group.Get(_.tracking_id/gslice, _.external_scaled_frags, _.sample_names, _.label_names, _.condition, _.replicate).Copy()

#second part of manual implementation of inverse HArray
names      = read_group.sample_names.Each(str.lower)()
dataraw    = [read_group[...,i].raw_frags/(','.join(names[i]) + '_raw') for i in range(len(names))]
datanorm   = [read_group[...,i].external_scaled_frags/(names[i] + '_normalized') for i in range(len(names))]
read_group = read_group.Get(_.tracking_id/gslice, *(dataraw + datanorm)).Copy()

#now get _exp.diff data for each comparison
results = [read_group];

all_expdiff = Read(files[diff_data_index]).Copy()
all_expdiff = all_expdiff / ('test_id', 'gene_id', 'gene', 'locus', 'sample_1', 'sample_2', 'status', 'value_1', 'value_2', 'log2_fold_change', 'test_stat', 'p_value', 'q_value', 'significant');
all_expdiff = all_expdiff.To(_.test_id,          Do=_.Cast(bytes));
all_expdiff = all_expdiff.To(_.gene_id,          Do=_.Cast(bytes));
all_expdiff = all_expdiff.To(_.gene,             Do=_.Cast(bytes));
all_expdiff = all_expdiff.To(_.locus,            Do=_.Cast(bytes));
all_expdiff = all_expdiff.To(_.sample_1,         Do=_.Cast(int));
all_expdiff = all_expdiff.To(_.sample_2,         Do=_.Cast(int));
all_expdiff = all_expdiff.To(_.status,           Do=_.Cast(bytes));
all_expdiff = all_expdiff.To(_.value_1,          Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.value_2,          Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.log2_fold_change, Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.test_stat,        Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.p_value,          Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.q_value,          Do=_.Cast(float));
all_expdiff = all_expdiff.To(_.significant,      Do=_.Cast(bytes));
all_expdiff = all_expdiff.Copy();

for i in xrange(len(C.cuffdiff_cmp)):
  a, b = C.cuffdiff_cmp[i];
  cname = (C.label_names[a] + "_" + C.label_names[b]).lower()
  expdiff = all_expdiff[(_.sample_1==a) & (_.sample_2==b)].Copy()

  expdiffd = expdiff.Without(_.gene_id, _.test_id, _.locus, _.gene, _.sample_1, _.sample_2)
  expdiffd = expdiffd / tuple([cname + '_' + name for name in expdiffd.Names])
  res = expdiff.Get(gslice, expdiffd)

  results.append(res.Copy());
#efor

#combine the data
R = results[0];
for i in range(1, len(results)):
    R = (R |Match(gslice, gslice)| results[i]).Copy()
#efor

#if annots file is available, add that one too (first slice should be id slice).
for annot_file in C.annotation_files:
    annots = Load(annot_file)
    R = annots |Match(gslice,gslice)| R

R = R % C.cuffdiff_test_type;

res2 = R |Match| read_group_compact
res2 = res2 % C.cuffdiff_test_type;

Save(res2, C.__cuffdiff_combine_output__()[0]);
Save(R, C.__cuffdiff_combine_output__()[1]);

sys.exit(0);
