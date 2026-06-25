



import matplotlib.pyplot as plt
import numpy as np
from . import const




def parse_poni(poni_str):
    raw_data = {}
    for line in poni_str.strip().split("\n"):
       
        if line.startswith("#") or not line:
            continue
        key, value = line.split(":", 1)
        raw_data[key.strip()] = value.strip()
        
    return  float(raw_data["Poni1"]), float(raw_data["Poni2"]), float(raw_data["Distance"]), float(raw_data["Wavelength"]),





def make_img_stats(self):
    self.img_mean = self.imgs.mean(axis=0)
    self.img_std = self.imgs.std(axis=0)
    
    self.img_mask = (self.imgs[0]!=4294967295)

    xmin =  - self.det_center[1] 
    xmax = self.det_shape[1]*self.det_px - self.det_center[1]
    ymin = - self.det_center[0]
    ymax = self.det_shape[0]*self.det_px - self.det_center[0]
    
    x_array = np.linspace(xmin,xmax, num=self.det_shape[1])
    y_array = np.linspace(ymin, ymax, num=self.det_shape[0])
    xx, yy = np.meshgrid(x_array, y_array)
    r = np.sqrt( xx**2 + yy**2)
    q = ((4*np.pi)/self.det_lam)*np.sin(0.5*np.arctan2(r, self.det_d))/1e10

    self.img_x = xx
    self.img_y = yy
    self.img_r = r
    self.img_q = q
    self.img_phi = np.arctan2(yy,xx)

    

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
    
    
    
    
    
    