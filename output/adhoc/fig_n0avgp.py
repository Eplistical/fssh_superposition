import numpy as np
from matplotlib import pyplot as plt
import sys
from templates import *

# parameters
savefname = 'fig_n0avgp.eps' 


# load and plot
fig, axes = plt.subplots(3,2, figsize=(10,8))
plt.subplots_adjust(top=0.9, bottom=0.18, left=0.1, right=0.9,hspace=0, wspace=0.25)

for i, A in enumerate(('0.02', '0.05', '0.10')):
    e_raw = np.loadtxt('e.A%s.scan_s.out' % (A))
    f = np.loadtxt('f.A%s.scan_s.out' % (A))
    fc = np.loadtxt('fc.A%s.scan_s.out' % (A))
    eh = np.loadtxt('eh.A%s.scan_s.out' % (A))

    # exact -- diab repr
    e = np.copy(e_raw)
    e[:,1], e[:,3] = e_raw[:,3], e_raw[:,1]
    e[:,2], e[:,4] = e_raw[:,4], e_raw[:,2]
    e[:,5], e[:,7] = e_raw[:,7], e_raw[:,5]
    e[:,6], e[:,8] = e_raw[:,8], e_raw[:,6]

    plot_n0(axes[i,0], e, f, fc, eh)
    plot_avgp(axes[i,1], e, f, fc, eh)

# annotation
xy1=(0.6,0.1)
axes[0,1].annotate('$A = 0.02$', xy=xy1, xycoords='axes fraction')
axes[1,1].annotate('$A = 0.05$', xy=xy1, xycoords='axes fraction')
axes[2,1].annotate('$A = 0.10$', xy=xy1, xycoords='axes fraction')

# y lim, tick and labels
set_y_label(axes[:,1], ['27', '30', '33'])
set_y_label(axes[:,0], ['0', '0.25', '0.5', '0.75', '1'])

# x lim, tick and labels
set_x_label(axes[:,0], ['0','0.25','0.5','0.75','1'])
set_x_label(axes[:,1], ['0','0.25','0.5','0.75','1'])


axes[2,1].set_xlabel('$c_1^2$')
axes[2,0].set_xlabel('$c_1^2$')

# title and legend
axes[0,1].set_title('Average Momentum')
axes[0,0].set_title('Adiabat 0 Population')

axes[2,0].legend(loc='upper left', bbox_to_anchor=(0.02, -0.3), ncol=4)

# save
fig.savefig(savefname)
