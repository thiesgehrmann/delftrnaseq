import numpy
import pylab as pl

def autolabel(rects, ax):
    # attach some text labels
    maxheight = -numpy.inf
    for rect in rects:
        height = rect.get_height()
        maxheight = max(height, maxheight)    
        if height < 1.0:
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%.2g'%height,
                    ha='center', va='bottom')
        else:
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                    ha='center', va='bottom')
    y = ax.get_ylim()
    if y[1] < maxheight * 1.1:
        ax.set_ylim([y[0], maxheight * 1.1])

def dual_bargraph(data1, data2, ylabel, xticklabels, legendlabels, figsize=None):
    if figsize is None:
        fig = pl.figure()
    else:
        fig = pl.figure(figsize=figsize)

    ax = fig.add_subplot(111)

    ind = numpy.arange(len(xticklabels)) + 1
    width = 0.35

    rects1 = ax.bar(ind-width, data1, width, color='r')
    rects2 = ax.bar(ind, data2, width, color='b')

    ax.set_ylabel(ylabel)
    ax.set_xticks(ind)
    ax.set_xticklabels(xticklabels, rotation=17)
    ax.legend( (rects1[0], rects2[0]), legendlabels,loc='best')

    autolabel(rects1, ax)
    autolabel(rects2, ax)
    fig.tight_layout()
    return fig


