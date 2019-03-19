import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys

plt.rcParams.update({'font.size': 15})

default_labels = [
r'exact adiabat 0',
r'exact adiabat 1',
r'FSSH adiabat 0',
r'FSSH adiabat 1',
r'FSSH phase-corr adiabat 0',
r'FSSH phase-corr adiabat 1',
]

def set_x_label(axes, ticklabels):
    ticks = [float(x) for x in ticklabels]
    lims = [float(ticklabels[0]), float(ticklabels[-1])]
    N = axes.size
    ticklabels2 = ['' for i in range(N)]

    for i, ax in enumerate(axes):
        ax.set_xlim(lims)
        ax.set_xticks(ticks)
        if i == N-1:
            ax.set_xticklabels(ticklabels)
        else:
            ax.set_xticklabels(ticklabels2)



def set_y_label(axes, ticklabels):
    ticklabels2 = [x for x in ticklabels]
    ticklabels2[0] = ''
    ticks = [float(x) for x in ticklabels]
    lims = [float(ticklabels[0]), float(ticklabels[-1])]
    N = axes.size

    for i, ax in enumerate(axes):
        ax.set_ylim(lims)
        ax.set_yticks(ticks)
        if i == N-1:
            ax.set_yticklabels(ticklabels)
        else:
            ax.set_yticklabels(ticklabels2)


def plot_trans_px_n(ax, e, f, fc, labels=default_labels):
    ax.plot(e[:,0], e[:,1], marker='s', color='green', linestyle='solid', linewidth=2, label=labels[0])
    ax.plot(e[:,0], e[:,3], marker='s', color='#4682B4', linestyle='solid', linewidth=2, label=labels[1])

    ax.plot(f[:,0], f[:,1], marker='', color='orange', linestyle='dashed', linewidth=2, label=labels[2])
    ax.plot(f[:,0], f[:,3], marker='', color='#BA55D3', linestyle='dashed', linewidth=2, label=labels[3])

    ax.plot(fc[:,0], fc[:,1], marker='', color='black', linestyle='dotted', linewidth=2, label=labels[4])
    ax.plot(fc[:,0], fc[:,3], marker='', color='red', linestyle='dotted', linewidth=2, label=labels[5])

    ax.tick_params(axis='both', direction='in', left=True, right=True, top=True, bottom=True)



def plot_trans_px_px(ax, e, f, fc, labels=default_labels):
    ax.plot(e[:,0], e[:,5], marker='s', color='green', linestyle='solid', linewidth=2, label=labels[0])
    ax.plot(e[:,0], e[:,7], marker='s', color='#4682B4', linestyle='solid', linewidth=2, label=labels[1])

    ax.plot(f[:,0], f[:,5], marker='', color='orange', linestyle='dashed', linewidth=2, label=labels[2])
    ax.plot(f[:,0], f[:,7], marker='', color='#BA55D3', linestyle='dashed', linewidth=2, label=labels[3])

    ax.plot(fc[:,0], fc[:,5], marker='', color='black', linestyle='dotted', linewidth=2, label=labels[4])
    ax.plot(fc[:,0], fc[:,7], marker='', color='red', linestyle='dotted', linewidth=2, label=labels[5])
