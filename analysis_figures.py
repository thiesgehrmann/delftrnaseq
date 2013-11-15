from ibidas import *
from sklearn.decomposition import PCA
from sklearn.cluster import MeanShift, estimate_bandwidth, Ward, KMeans
from sklearn import metrics
import pylab as pl
import numpy
from collections import Counter
from matplotlib_venn import venn2, venn3
from figure_tools import *

markers = ['o', '1', '8','2', 'p', '3', 's', '4', 'h','H','D','d','*','+', 'x', '.']


def create_pca(data, filename, logfilter=False):
    x = data.external_scaled_frags()
    #print x.shape
    #print numpy.std(numpy.log2(x+1), axis=1)
    #x = x[numpy.std(numpy.log2(x+50),axis=1) > 1.0,:]
    #print x.shape

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
    names = Rep(data.Names)[_.HasPattern('significant')]()
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

    fig = dual_bargraph(up, down, "Number of significant genes", cnames, \
                  ('upregulated','downregulated'),figsize=(15.0, 9.0))
    pl.savefig(filename, dpi=200)


def create_venn(data, compare_sets, filename):
    compare_names = []
    for compare_set in compare_sets:
        label_names = [data[_.condition == condid].label_names[0]() for condid in compare_set]
        compare_names.append(("_".join(label_names)).lower())

    indicators = []
    for compare_name in compare_names:
        indicators.append((data.Get(compare_name + "_significant") == "yes")())


    skipvenn = any([numpy.sum(filter)==0 for filter in indicators])
    
    if len(compare_names) == 2:
        indicators = ["%d%d" % x for x in zip(*indicators)]
    else:
        assert len(compare_names) == 3, 'Number of compare sets in Venn diagram should be 2 or 3'
        indicators = ["%d%d%d" % x for x in zip(*indicators)]

    subsets = Counter(indicators)
    
    fig = pl.figure()
    fig.set_facecolor('white')
    ax = fig.add_subplot(111)
    
    if not skipvenn:
        if len(compare_names) == 2:
            venn2(subsets, set_labels=tuple(compare_names),ax=ax)
        else:
            venn3(subsets, set_labels=tuple(compare_names),ax=ax)
    else:
        pl.text(0.5, 0.5, 'One or more sets has no significant genes', horizontalalignment='center')
        pl.axis('off')

    pl.savefig(filename, dpi=200)
   



