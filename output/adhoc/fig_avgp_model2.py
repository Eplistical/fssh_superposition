import numpy as np
from matplotlib import pyplot as plt
import sys
from templates import *

# parameters
savefname = 'fig_avgp_model2.eps' 


# load and plot
fig, axes = plt.subplots(3,1, figsize=(6.5,8))
plt.subplots_adjust(top=0.9, bottom=0.18, left=0.1, right=0.9,hspace=0, wspace=0.25)

for i, A in enumerate(('0.02', '0.05', '0.10')):
    e_raw = np.loadtxt('e.A%s.scan_s.model2.out' % (A))
    f = np.loadtxt('f.A%s.scan_s.model2.out' % (A))
    fc = np.loadtxt('fc.A%s.scan_s.model2.out' % (A))
    eh = np.loadtxt('eh.A%s.scan_s.model2.out' % (A))

    # exact -- diab repr
    e = np.copy(e_raw)
    e[:,1], e[:,3] = e_raw[:,3], e_raw[:,1]
    e[:,2], e[:,4] = e_raw[:,4], e_raw[:,2]
    e[:,5], e[:,7] = e_raw[:,7], e_raw[:,5]
    e[:,6], e[:,8] = e_raw[:,8], e_raw[:,6]

    plot_avgp(axes[i], e, f, fc, eh)

# annotation
xy1=(0.6,0.1)
axes[0].annotate('$A = 0.02$', xy=xy1, xycoords='axes fraction')
axes[1].annotate('$A = 0.05$', xy=xy1, xycoords='axes fraction')
axes[2].annotate('$A = 0.10$', xy=xy1, xycoords='axes fraction')

# y lim, tick and labels
set_y_label(axes, ['27', '30', '33'])

# x lim, tick and labels
set_x_label(axes, ['0','0.25','0.5','0.75','1'])


axes[2].set_xlabel('$c_1^2$')

# title and legend
axes[0].set_title('Average Momentum')

axes[2].legend(loc='upper left', bbox_to_anchor=(0.05, -0.29), ncol=2)

# save
fig.savefig(savefname)
