import matplotlib.pyplot as plt
import numpy as np
from . import const



    

def plot_iq(self, yerrs=None, xlim=None, ylim=None, logX=None, logY=None, **kwargs):
    if self.det=='eiger':
        if xlim is None:
            xlim = const.EIGER_QRANGE
        if logX is None: 
            logX = True
        if logY is None: 
            logY = True
        
    if self.det=='pilatus':
        if xlim is None:
            xlim = const.PILATUS_QRANGE
        if logX is None: 
            logX = False
        if logY is None: 
            logY = True
            
    if yerrs is not None:
        yerrs = self.Is_err
    
    _plot_1d(self.qs, self.Is, yerrs=yerrs, xlim=xlim, ylim=ylim, logX=logX, logY=logY, **kwargs)




def _plot_1d(x, ys, yerrs=None, new_fig=True, logX=None, logY=None, cmap='viridis', xlim=None, ylim=None, mask=None):
    
    if new_fig:
        plt.figure()

    if yerrs is None:
        yerrs = np.zeros(ys.shape)

    if mask is None:
        mask = np.array([True]*ys.shape[0])

    if type(cmap) is str:
        cmap = plt.cm.get_cmap(cmap, ys.shape[0])
    
        
    for i, (y, yerr, mask_val) in enumerate(zip(ys, yerrs, mask)):
        if mask_val:
            plt.errorbar(x, y, yerr=yerr, color=cmap(i), label=i)
    plt.xlabel('q [1/A]')
    plt.ylabel('Inten.')
    
    if logX: plt.xscale('log')
    if logY: plt.yscale('log')
    if xlim is not None: plt.xlim(xlim)
    if ylim is not None: plt.ylim(ylim)
        

def _plot_2d(img, log=False, vmin=None, vmax=None, det_center=(0,0), det_px = 1):

    if log:
        img=np.log10(np.abs(img)+1)
        

    xmin =  - det_center[1] 
    xmax = img.shape[1]*det_px - det_center[1]
    ymin = - det_center[0]
    ymax = img.shape[0]*det_px - det_center[0]
    extent_centered = [xmin, xmax, ymin, ymax]
    
    plt.figure()
    plt.imshow(img, extent=extent_centered, cmap='viridis', clim=(vmin,vmax), origin='lower')
    plt.colorbar()
    plt.scatter(0,0,marker='+', color='red')
    
    
    