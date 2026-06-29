import matplotlib.pyplot as plt
import numpy as np
from . import const



    

def plot_iq(self, i=0, axes=None, color='red', yerr=None, xlim=None,logX=True, logY=True, **kwargs):
    
    
    if axes is None:
        plt.figure()
        axes = plt.gca()

    if yerr is None:
        yerr = np.zeros(self.qs.shape)

    
    
    axes.errorbar(self.qs, self.Is[i], yerr=yerr, color=color,**kwargs)
          
    axes.set_xlabel('q [1/A]')
    axes.set_ylabel('Inten.')
    
    if logX: axes.set_xscale('log')
    if logY: axes.set_yscale('log')

    if xlim is None:
        xlim = self.det.qrange
    axes.set_xlim(xlim)
   



def plot_iqs(self, axes=None, cmap='viridis', xlim=None, ylim=None, logX=True, logY=True, color=None, **kwargs):

    if axes is None:
        plt.figure()
        axes = plt.gca()



    if cmap is not None:
        cmap = plt.colormaps.get_cmap(cmap)
    else:
        def cmap(x):
            return color

    for i in range(self.n_ims):
        axes.plot(self.qs, self.Is[i], color=cmap(i / (self.n_ims - 1)), **kwargs)

    axes.set_xlabel('q [1/A]')
    axes.set_ylabel('Inten.')
    
    if logX: axes.set_xscale('log')
    if logY: axes.set_yscale('log')

    if xlim is None:
        xlim = self.det.qrange
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

    
    
    
    
    
    
    



    
    