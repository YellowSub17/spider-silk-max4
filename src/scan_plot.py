

import matplotlib.pyplot as plt
import numpy as np
from .const import EIGER_QRANGE, PILATUS_QRANGE

class ScanPlot:


    def plot_iq(self, errs=False, logX=None,logY=None, xlim=None):            
        self._plot(self.qs, self.Is, logX=logX,logY=logY,xlim=xlim)
       

    def plot_iq_errs(self,logX=None,logY=None, xlim=None):
        self._plot(self.qs, self.Is, self.Is_err, logX=logX,logY=logY,xlim=xlim)

    

    def _plot(self, x, ys, yerrs=None, new_fig=True, logX=None, logY=None, cmap='viridis', xlim=None):
        
        if new_fig:
            plt.figure()
    
        if yerrs is None:
            yerrs = np.zeros(ys.shape)
            
        
                
        

        cmap = plt.cm.get_cmap('viridis', ys.shape[0])
        for i, (y, yerr) in enumerate(zip(ys, yerrs)):
            plt.errorbar(x, y, yerr=yerr, color=cmap(i))

        plt.xlabel('q [1/A]')
        plt.ylabel('Inten.')
        if self.det=='eiger':
            if xlim is None:
                plt.xlim(EIGER_QRANGE)
            if logX is None: 
                logX = True
            if logY is None: 
                logY = True
                
        if self.det=='pilatus':
            if xlim is None:
                plt.xlim(PILATUS_QRANGE)
            if logX is None: 
                logX = False
            if logY is None: 
                logY = True
        
        
        if logX: plt.xscale('log')
        if logY: plt.yscale('log')
            
    
    
    def plot_img(self, i=0, log=False, xlim=None, ylim=None):
        
        if i is None:
            img = np.mean(self.get_imgs(), axis=0)
        else:
            img = self.get_imgs()[i]

        if log:
            img=np.log10(img)
            
        extent_centered = [
            0 - self.det_center[1],          # xmin
            img.shape[1]*self.det_px - self.det_center[1],    # xmax
            img.shape[0]*self.det_px - self.det_center[0],   # ymin (bottom)
            0 - self.det_center[0],           # ymax (top),
        ]
        
        plt.figure()
        plt.imshow(img, extent=extent_centered,cmap='viridis')
        plt.colorbar()
        plt.scatter(0,0,marker='+', color='red')
        if xlim is None:
            xlim = [0 - self.det_center[1],img.shape[1]*self.det_px - self.det_center[1]]
            plt.xlim(xlim)
        if ylim is None:
            ylim = [0 - self.det_center[0],img.shape[0]*self.det_px - self.det_center[0]]
            plt.ylim(ylim)