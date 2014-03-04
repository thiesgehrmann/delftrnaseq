import numpy
import pylab as pl

def autolabel(rects, ax) :
    # attach some text labels
    maxheight = -numpy.inf
    for rect in rects:
        height = rect.get_height()
        maxheight = max(height, maxheight)    
        #if height < 1.0:
        if isinstance(height, float) and height < 1.0 :
            ax.text(rect.get_x() + rect.get_width() / 2.0 , 1.05 * height, '%.2g' % height, ha = 'center', va = 'bottom', fontsize = 18)
        elif isinstance(height, float) :
            ax.text(rect.get_x() + rect.get_width()/2., 1.05 * height, '%.2f' % height, ha = 'center', va = 'bottom', rotation = 90, fontsize = 18)
        else:
            ax.text(rect.get_x() + rect.get_width()/2., 1.05 * height, format(int(round(height)), ",d"), ha = 'center', va = 'bottom', rotation = 90, fontsize = 18)
    y = ax.get_ylim()
    if y[1] < maxheight * 1.2 :
        ax.set_ylim([y[0], maxheight * 1.3])

def dual_bargraph(data1, data2, ylabel, xticklabels, legendlabels, figsize = None) :
    if figsize is None :
        fig = pl.figure()
    else :
        fig = pl.figure(figsize = figsize)

    ax = fig.add_subplot(111)

    ind = numpy.arange(len(xticklabels)) + 1
    width = 0.65
    scaler = 2.0

    rects1 = ax.bar(scaler * ind - width, data1, width, color = 'r')
    rects2 = ax.bar(scaler * ind, data2, width, color = 'b')

    x = ax.get_xlim()
    ax.set_xlim([x[0], numpy.max(ind) * scaler + width + (numpy.min(ind) * scaler - width - x[0])])
    ax.set_ylabel(ylabel, fontsize = 20)
    ax.set_xticks(scaler * ind)
    ax.set_xticklabels(xticklabels, rotation = 17)
    ax.legend( (rects1[0], rects2[0]), legendlabels, loc = 'best', fontsize = 18)

    autolabel(rects1, ax)
    autolabel(rects2, ax)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 18)
    fig.tight_layout()
    return fig

def tripple_bargraph(data1, data2, data3, ylabel, xticklabels, legendlabels, figsize = None):
    if figsize is None :
        fig = pl.figure()
    else :
        fig = pl.figure(figsize = figsize)

    ax = fig.add_subplot(111)

    ind = numpy.arange(len(xticklabels)) + 1
    width = 0.65
    scaler = 3.0

    rects1 = ax.bar(scaler * ind - width, data1, width, color = 'r')
    rects2 = ax.bar(scaler * ind, data2, width, color = 'b')
    rects3 = ax.bar(scaler * ind + width, data3, width, color = 'y')

    x = ax.get_xlim()
    ax.set_xlim([x[0], numpy.max(ind) * scaler + 2 * width + (numpy.min(ind) * scaler - width - x[0])])
    ax.set_ylabel(ylabel, fontsize = 20)
    ax.set_xticks(scaler * ind + width / 2.0)
    ax.set_xticklabels(xticklabels, rotation = 17)
    ax.legend( (rects1[0], rects2[0], rects3[0]), legendlabels, loc = 'best', fontsize = 18)

    autolabel(rects1, ax)
    autolabel(rects2, ax)
    autolabel(rects3, ax)
    ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 18)
    fig.tight_layout()
    return fig
