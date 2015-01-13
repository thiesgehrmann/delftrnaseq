from ibidas import *
from sklearn.decomposition import PCA
from sklearn.cluster import MeanShift, estimate_bandwidth, Ward, KMeans
from sklearn import metrics
import pylab as pl
import numpy
from collections import Counter
from figure_tools import *

from matplotlib import pyplot as plt;
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles;
  
markers = ['o', '1', '8','2', 'p', '3', 's', '4', 'h','H','D','d','*','+', 'x', '.']


def create_pca(data, filename, logfilter=False):
    x = data.external_scaled_frags()

    condlabels = data.label_names()
    samplelabels = data.sample_names()

    if logfilter:
        x = numpy.log2(x + 16)
        x = x[numpy.std(x,axis=1) > 0.5,:]

    pca = PCA(n_components=2)
    X_r = pca.fit(x.T).transform(x.T)
    ratio = pca.explained_variance_ratio_

    fig = pl.figure()
    ax = pl.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ucondlabels = numpy.unique(condlabels)
    for i, clabel in enumerate(ucondlabels):
        t = X_r[condlabels == clabel,:]
        pl.scatter(t[:,0], t[:,1], c=pl.cm.jet(float(i)/(len(ucondlabels)-1)), s=250, label=clabel, alpha=0.7)

    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    pl.title('PCA plot of all samples')
    pl.xlabel('Explained variance: {0:.2f}'.format(ratio[0]))
    pl.ylabel('Explained variance: {0:.2f}'.format(ratio[1]))

    pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    pl.savefig(filename,dpi=200)
    #pl.show() 

def create_pca2(data, filename, logfilter=False):
    labels  = dict([ (n.split('_')[int(n.split('_')[-1])-1], n) for n in data.Names if '_value_' in n ]).items()
    names   = [ i[1] for i in labels ];
    samples = [ i[0] for i in labels ];

    x = numpy.array(zip( *data.Get(*names)() ));
    
    if logfilter:
        x = numpy.log2(x + 16)
        x = x[numpy.std(x,axis=1) > 0.5,:]
    
    pca = PCA(n_components=2)
    X_r = pca.fit(x.T).transform(x.T)
    ratio = pca.explained_variance_ratio_
    
    fig = pl.figure()
    ax = pl.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    for i, clabel in enumerate(samples):
        t = X_r[i,:]
        pl.scatter(t[0], t[1], c=pl.cm.jet(float(i)/(len(samples)-1)), s=250, label=clabel, alpha=0.7)
    
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    pl.title('PCA plot of all samples')
    pl.xlabel('Explained variance: {0:.2f}'.format(ratio[0]))
    pl.ylabel('Explained variance: {0:.2f}'.format(ratio[1]))
    
    pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    pl.savefig(filename,dpi=200)
    #pl.show() 


def create_cluster(data, filename):
    x = data.external_scaled_frags()
    x = numpy.log2(x + 16)

    #bandwidth = estimate_bandwidth(x, quantile=0.2, n_samples=500)
    #ms = MeanShift(bandwidth = bandwidth, bin_seeding=True)
    ms = KMeans(n_clusters=12)

    x = x[numpy.std(x, axis=1) > 0.5,:]
    x = (x - numpy.mean(x, axis=1)[:, numpy.newaxis])
    x = x / (numpy.std(x,axis=1)[:,numpy.newaxis] + 1e-10)


    ms.fit(x)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_

    nclusters = len(cluster_centers)

    rows = numpy.ceil(numpy.sqrt(nclusters))
    cols = numpy.ceil(nclusters / rows)
   
    fig = pl.figure(figsize=(20, 12.0))

    xticks = numpy.arange(x.shape[1])
    for k, col in enumerate(range(nclusters)):
        class_members = labels == k
        pl.subplot(rows, cols, k+1)
        ty = x[class_members,:]
        tx = numpy.tile(xticks, ty.shape[0]).reshape(ty.shape[0],len(xticks))
        pl.plot(tx.T, ty.T,'k',alpha=0.05)
        pl.plot(xticks, cluster_centers[k,:],'r',linewidth=3)
        pl.xlim([0, len(xticks)-1])
        pl.title('Cluster ' + str(k+1) + ' (%d genes)' % ty.shape[0])
        if ((k+1) / cols) > 3:
            pl.xticks(xticks, data.sample_names(), rotation=30, ha='right')
        pl.ylabel('Normalized Fold Change')
    pl.savefig(filename, dpi=200)        
    #pl.show()


def create_diffgenes_stats(data, filename):
    names = [n for n in data.Names if "_significant" in n ]
    up = []
    down = []
    cnames = []
    for name in names:
        comparison_name = name[:-12]
        xdata = data[_.Get(name) == "yes"]
        countup = (xdata.Get(comparison_name + "_log2_fold_change") > 0).Sum()()
        countdown = (xdata.Get(comparison_name + "_log2_fold_change") < 0).Sum()()
        cnames.append(comparison_name)
        up.append(countup)
        down.append(countdown)

    print up, down
    fig = dual_bargraph(up, down, "Number of significant genes", cnames, \
                  ('upregulated','downregulated'),figsize=(15.0, 9.0))
    pl.savefig(filename, dpi=200)


#def create_venn(data, compare_sets, names, filename):
#    compare_names = []
#    for compare_set in compare_sets:
#        label_names = [ names[condid] for condid in compare_set]
#        compare_names.append(("_".join(label_names)).lower())
#
#    indicators = []
#    for compare_name in compare_names:
#        indicators.append((data.Get(compare_name + "_significant") == "yes")())
#
#    skipvenn = any([numpy.sum(filter)==0 for filter in indicators])
#    
#    if len(compare_names) == 2:
#        indicators = ["%d%d" % x for x in zip(*indicators)]
#    else:
#        assert len(compare_names) == 3, 'Number of compare sets in Venn diagram should be 2 or 3'
#        indicators = ["%d%d%d" % x for x in zip(*indicators)]
#
#    subsets = Counter(indicators)
#    
#    fig = pl.figure()
#    fig.set_facecolor('white')
#    ax = fig.add_subplot(111)
#    
#    if not skipvenn:
#        if len(compare_names) == 2:
#            venn2(subsets, set_labels=tuple(compare_names),ax=ax)
#        else:
#            venn3(subsets, set_labels=tuple(compare_names),ax=ax)
#    else:
#        pl.text(0.5, 0.5, 'One or more sets has no significant genes', horizontalalignment='center')
#        pl.axis('off')
#
#    pl.savefig(filename, dpi=200)

def create_venn(data, compare_sets, names, filenames, udsplit):

  datas = [];
  dirs  = [];

  compare_names = []
  for compare_set in compare_sets:
    label_names = [ names[condid] for condid in compare_set]
    compare_names.append(("_".join(label_names)).lower())
  #efor

  x = set(data.Names)
  if not 'gene_id' in x:
      name = 'test_id'
  else:
      name = 'gene_id'

  if udsplit:
    ups   = [];
    downs = [];
    for compare_name in compare_names:
      up   = set(data[data.Get((compare_name + '_significant')) == 'yes'].Get(name)()) & set(data[data.Get((compare_name + '_log2_fold_change')) > 0].Get(name)());
      down = set(data[data.Get((compare_name + '_significant')) == 'yes'].Get(name)()) & set(data[data.Get((compare_name + '_log2_fold_change')) < 0].Get(name)());
      ups.append(up);
      downs.append(down);
    #efor
    datas.append(ups);
    datas.append(downs);
    dirs = [ 'up', 'down' ];
  else:
    alls = [];
    for compare_name in compare_names:
      all = set(data[data.Get(compare_name + '_significant') == 'yes'].Get(name)());
      alls.append(all);
    #efor
    datas.append(alls);
    dirs = [ 'all' ];
  #fi

  numbers  = [];
  overlaps = [];
  for (data, filename, dir) in zip(datas, filenames, dirs):
    if len(data) == 2:
      numbers = [0]*3;
      titles  = [ compare_names[0], compare_names[1], '%s n %s' % (compare_names[0], compare_names[1]) ];
      numbers[2] = data[0] & data[1];
      numbers[0] = data[0] - numbers[2];
      numbers[1] = data[1] - numbers[2]; 
    elif len(data) == 3:
      numbers = [0]*7;
      titles = [ compare_names[0], compare_names[1], ' n '.join([compare_names[0], compare_names[1]]), compare_names[2], ' n '.join([compare_names[0], compare_names[2]]), ' n '.join([compare_names[1], compare_names[2]]), ' n '.join([compare_names[0], compare_names[1], compare_names[2]]) ];
      numbers[6] = data[0] & data[1] & data[2];
      numbers[5] = data[1] & data[2] - numbers[6];
      numbers[4] = data[0] & data[2] - numbers[6];
      numbers[3] = data[2] - numbers[4] - numbers[5] - numbers[6];
      numbers[2] = data[0] & data[1] - numbers[6];
      numbers[1] = data[1] - numbers[2] - numbers[5] - numbers[6];
      numbers[0] = data[0] - numbers[2] - numbers[4] - numbers[6];
    else:
      print "Error.";
    #fi
    overlaps.extend([ (','.join(compare_names) + ' ' + dir, t, list(n)) for (t,n) in zip(titles,numbers)]);
    draw_venn(','.join(compare_names) + ' ' + dir, compare_names, [len(s) for s in numbers], filename);
  #efor

  return overlaps;
#edef

###############################################################################

def draw_venn(title, names, numbers, out):
  
  if len(numbers) == 7:
    if numbers[0] + numbers[2] + numbers[4] + numbers[6] == 0:
      numbers = [ numbers[1], numbers[3], numbers[5] ];
      names   = [ names[1], names[2] ];
    elif numbers[1] + numbers[2] + numbers[5] + numbers[6] == 0:
      numbers = [ numbers[0], numbers[3], numbers[4] ];
      names   = [ names[0], names[2] ];
    elif numbers[3] + numbers[4] + numbers[5] + numbers[6] == 0:
      numbers = [ numbers[0], numbers[1], numbers[2] ];
      names   = [ names[0], names[1] ];
    #fi
  #fi
  
  plt.cla();
  plt.figure(figsize=(10,10))
  if len(numbers) == 7:
    plt.cla();
    plt.figure(figsize=(10,10))
    v = venn3(subsets=numbers, set_labels = names)
    c = venn3_circles(subsets=numbers, linestyle='dashed')
  else:
    v = venn2(subsets = numbers, set_labels = names);
    c = venn2_circles(subsets = numbers, linestyle='dashed');
  #fi
  
  plt.title(title)
  plt.savefig(out);
#edef
