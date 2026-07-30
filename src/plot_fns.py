import matplotlib.pyplot as plt
import numpy as np
from . import const
from . import det

###
# Plotting functions that are shared between the scan and comboscan classes
###


def plot_1D(self, x,y, xlabel='q [1/A]', ylabel='Inten.', axes=None, color='red', yerr=None, xlim=None,logX=True, logY=True, **kwargs):
    if axes is None:
        plt.figure()
        axes = plt.gca()

    
    if yerr is not None and np.any(yerr):
        axes.errorbar(x, y, yerr=yerr, color=color, **kwargs)
    else:
        axes.plot(x, y, color=color, **kwargs)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    if logX: axes.set_xscale('log')
    if logY: axes.set_yscale('log')

    if xlim is not None:
        axes.set_xlim(xlim)


def plot_1Ds(self, x, ys, xlabel='q [1/A]', ylabel='Inten.', axes=None, cmap='viridis', xlim=None, ylim=None, logX=True, logY=True, color=None, **kwargs):

    '''
    Plots multiple merged saxs and waxs profiles
    '''

    if axes is None:
        plt.figure()
        axes = plt.gca()

    if color is not None:
        def cmap(x):
            return color
    else:
        cmap = plt.colormaps.get_cmap(cmap)

    for i in range(self.n_ims):
        axes.plot(x, ys[i], color=cmap(i / (self.n_ims - 1)), **kwargs)

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)

    if logX: axes.set_xscale('log')
    if logY: axes.set_yscale('log')

    if xlim is not None:
        axes.set_xlim(xlim)


    





def plot_img(self, i=0, log=False, vmin=None, vmax=None, mask=True):
    img = self.imgs[i].astype(float)
    print(img.dtype)
    if log:
        img = np.log10(np.abs(img)+1)
    if mask:
        img[~self.det.mask] = np.nan

    extent_centered = [self.det.xrange[0], self.det.xrange[1], self.det.yrange[0], self.det.yrange[1]]
    plt.figure()
    plt.imshow(img, extent=extent_centered, cmap='viridis', clim=(vmin,vmax), origin='lower')
    plt.colorbar()
    plt.scatter(0,0,marker='+', color='red')













