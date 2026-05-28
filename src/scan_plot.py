

import matplotlib.pyplot as plt
import numpy as np

class ScanPlot:


    def plot_iq(self, logI=False, logq=False, errs=False):
        ## todo: add option to plot either or both the detectors
        self._plot(self.qs, self.Is)
        
        #if errs:
        #    plot_fn = plt.errorbars
        #else:
        #    plot_fn = plt.plot

        #cmap = plt.cm.get_cmap('viridis', self.n_ims)
        #plt.figure()
        #q = self.qs[:]
        #if logq:
        #    q = np.log10(q)
        #for i in range(self.n_ims):
        #    I = self.Is[i,:]
        #    if logI:
        #        I=np.log10(I)
        #    plt.plot(q, I, color=cmap(i))
        #    
        #plt.xlabel('q [1/A]')
        #plt.ylabel('Inten.')

    def plot_iq_errs(self):
        ## todo: add option to plot either or both the detectors

        self._plot(self.qs, self.Is, self.Is_err)
        #cmap = plt.cm.get_cmap('viridis', self.n_ims)
        ##plt.figure()
        #for i in range(self.n_ims):
        #    plt.errorbar(self.qs, self.Is[i,:], self.Is_err[i,:], capsize=1, linestyle='', color=cmap(i))

        



    def _plot(self, x, ys, yerrs=None, new_fig=True, logX=False, logY=True, cmap='viridis'):
        if new_fig: plt.figure()

        if yerrs is None:
            yerrs = np.zeros(ys.shape)

        cmap = plt.cm.get_cmap('viridis', ys.shape[0])
        for i, (y, yerr) in enumerate(zip(ys, yerrs)):
            plt.errorbar(x, y, yerr=yerr, color=cmap(i))

        
        
        plt.xlabel('q [1/A]')
        plt.ylabel('Inten.')
        
        if logX: plt.xscale('log')
        if logY: plt.yscale('log')
            
        



    #plt.errorbar(q,I, yerr=errors, capsize=1, linestyle='')
    
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