import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import matplotlib
import pylab
from matplotlib.transforms import Affine2D
from matplotlib.collections import PathCollection
import matplotlib.pyplot as plt
import argparse
import sys

FULL_LABELS = ['Frontal', 'Parietal', 'Temporal', 'Occipital', 'Subcortical',
               'Occipital', 'Temporal', 'Parietal', 'Frontal']
TERR_LABELS = ['ACA','MCA','PCA','IT','PCA','MCA','ACA']
FULL_LABELS_IT = ['Frontal', 'Parietal', 'Temporal', 'Occipital',
                  'Subcortical/Infratentorial', 'Occipital', 'Temporal',
                  'Parietal', 'Frontal']
LABELS_LR = ['Frontal', 'Temporal', 'Subcortical', 'Occipital', 'Parietal']

ORDER=['FL_Layers1','PL_Layers1','TL_Layers1','OL_Layers1','BGIT_Layers1','OR_Layers1','TR_Layers1','PR_Layers1','FR_Layers1',
       'FL_Layers2','PL_Layers2','TL_Layers2','OL_Layers2','BGIT_Layers2','OR_Layers2','TR_Layers2','PR_Layers2','FR_Layers2',
       'FL_Layers3','PL_Layers3','TL_Layers3','OL_Layers3','BGIT_Layers3','OR_Layers3','TR_Layers3','PR_Layers3','FR_Layers3',
       'FL_Layers4','PL_Layers4','TL_Layers4','OL_Layers4','BGIT_Layers4','OR_Layers4','TR_Layers4','PR_Layers4','FR_Layers4']

def create_bullseye_plot(data, color, num_layers=4, num_lobes=9, vmin=0,
                         vmax=1, labels=FULL_LABELS, thr=None):
    n_layer_steps = num_layers + 1
    n_layer_steps_comp = n_layer_steps * 1j
    theta, rad = np.mgrid[0:2*np.pi:(num_lobes*40+1)*1j,
                          0.2:1:n_layer_steps_comp]
    print(theta.shape)
    z_value = np.tile(data, [40, 1]).T

    z_value = z_value.reshape([num_layers, num_lobes*40])

    z_new = np.asarray(np.concatenate((z_value, -10*np.ones([1, num_lobes*40])),
                                      axis=0)).T
    z_new = np.asarray(np.concatenate((z_new, -10*np.ones([1, n_layer_steps])),
                       axis=0))
    z_new = z_new.reshape(theta.shape)

    print(z_new.shape)
    rot = Affine2D().rotate_deg(30)

    colormap = plt.get_cmap(color)
    colormap.set_bad('white')
    if thr is not None:
        z_new = np.where(z_new < thr, -10*np.ones_like(z_new), z_new)
    z_new = np.ma.masked_values(z_new, -10)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    fig, axis = plt.subplots(ncols=1, subplot_kw=dict(projection='polar',
                                                      transform=rot))
    axis.pcolormesh(theta, rad, z_new, clip_on=True, cmap=colormap,alpha=1,
                    edgecolors='face', antialiased=True, linewidth=0,rasterized=True, vmin=vmin, vmax=vmax)
    for i in range(0, n_layer_steps):
        axis.plot(theta.T[0], np.repeat(rad[0][i]+0.1, theta.T[0].size), '-',
                color=[0.5, 0.5, 0.5], lw=1)
    for i in range(num_lobes):
        theta_i = i * 360/num_lobes * np.pi / 180
        axis.plot([theta_i, theta_i], [rad[0][0]-0.1, 0.9], '-',
                  color=[0.5, 0.5, 0.5], lw=1)

    axis.set_theta_zero_location("N")
    x_tab = np.arange(180/num_lobes*np.pi/180, 2*np.pi,
                      step=360/num_lobes*np.pi/180)
    x_lab = labels
    xdata = x_tab
    ydata = np.asarray([1, ]*len(x_tab))

    for x_value, y_value, l_value in zip(xdata, ydata, x_lab):
        if x_value > np.pi/2 and x_value < 3*np.pi/2:
            x_rot = x_value+np.pi
        else:
            x_rot = x_value
        pylab.text(x_value, y_value, l_value, rotation=np.rad2deg(x_rot),
                   fontsize=10,
                   horizontalalignment='center', verticalalignment='center')
    # for (t,l) in zip(xT,xL):
    #     plt.xticks(t, l, 0y=0.11, rotation=t)
    pylab.text(0.95,1.03,"Right",fontsize=10, transform=axis.transAxes)
    pylab.text(0.05, 1.03, "Left", fontsize=10, transform=axis.transAxes)
    axl = fig.add_axes([0.87, 0.1, 0.03, 0.8])
    cb1 = matplotlib.colorbar.ColorbarBase(axl, cmap=colormap, norm=norm,
                                           orientation='vertical')

    axis.set_yticklabels([])
    axis.set_xticklabels([])
    axis.grid(False)
    return plt.gcf()

def main(argv):

    parser = argparse.ArgumentParser(description='Assess number of lesions in file across files')
    parser.add_argument('-csv', dest='csv',required=True, metavar='input csv')
    parser.add_argument('-color', dest='color', metavar='colormap',
                        type=str, default='Reds',
                        help='color map to use')
    parser.add_argument('-type',default='FreqLobes',dest='type')
    parser.add_argument('-vmax',dest='vmax',default=0.3)
    
 


    try:
        args = parser.parse_args()
        # print(args.accumulate(args.integers))
    except argparse.ArgumentTypeError:
        print('print_be.py -csv XXX -color XXX ')

    df_les = pd.read_csv(args.csv)
    list_fields = [args.type+k for k in ORDER]
    data = np.asarray(df_les[list_fields])
    be = create_bullseye_plot(data,args.color,vmax=args.vmax)
    be.show()
    be.savefig(args.csv.split('.csv')[0]+'.png')

if __name__ == "__main__":
     args_temp = ['-csv','/Users/carolesudre/Data/MNI/test.csv',]
     #main(args_temp)
     main(sys.argv[1:])
