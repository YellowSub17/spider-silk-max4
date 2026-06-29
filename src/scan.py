
import numpy as np
import matplotlib.pyplot as plt
import h5py
import hdf5plugin

#from . import scan_plot 

from . import const

from . import plot_fns
from . import scale_fns

from .det import pilatus_det, eiger_det

class Scan:

    def __init__(self, scan_id=102470, det=pilatus_det, load_imgs=True, load_n=1):
        

        self.scan_id = scan_id
        self.det = det

        self.raw_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}.h5"
        self.raw_det_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}_{self.det.name}.h5"
        self.raw_pcap_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}_pcap.h5"
        self.azint_path = f"{const.DATA_PATH}/process/azint/scan-{self.scan_id}_{self.det.name}_integrated.h5"

        
        with h5py.File(self.azint_path,'r') as f:
            
            self.qs = f['/entry/data/radial_axis'][:]
            self.dq = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            
            if load_n=='all':
                self.Is = f['/entry/data/I'][...]
                self.Is_err = f['/entry/data/I_errors'][...]
            else:
                self.Is = f['/entry/data/I'][0:load_n,...]
                self.Is_err = f['/entry/data/I_errors'][0:load_n,...]

        
        
        with h5py.File(self.raw_path,'r') as f:
            if load_n=='all':
                self.n_ims = f[f'/entry/instrument/ring_current/data/'][...].shape[0]
                self.ring_current = f[f'/entry/instrument/ring_current/data'][...]
            else:
                self.n_ims = load_n
                self.ring_current = f[f'/entry/instrument/ring_current/data'][0:load_n,...]
            

        with h5py.File(self.raw_pcap_path,'r') as f:
            if load_n=='all':
                self.i0 = f['/entry/instrument/pandabox/data/i_0'][...]
                self.it = f['/entry/instrument/pandabox/data/i_t'][...]
            else:
                self.i0 = f['/entry/instrument/pandabox/data/i_0'][0:load_n,...]
                self.it = f['/entry/instrument/pandabox/data/i_t'][0:load_n,...]
                
        self.flags = (self.i0==self.i0)

        self.imgs = None
        if load_imgs:
            self.load_imgs()


    def load_imgs(self):
        if self.imgs is not None:
            return
        with h5py.File(self.raw_path,'r') as f:
            self.imgs = f[f'/entry/instrument/{self.det.name}/data/'][0:self.n_ims,...]
        self.img_mean = self.imgs.mean(axis=0)
        self.img_std = self.imgs.std(axis=0)
    
    def plot_img(self,i=0, **kwargs):
        return plot_fns.plot_img(self,i=i,**kwargs)

    def plot_iq(self, i=0, **kwargs):
        return plot_fns.plot_iq(self, i=i,**kwargs)

    def plot_iqs(self,**kwargs):
        return plot_fns.plot_iqs(self, **kwargs)


    
    def norm_qrange(self, qmin=0, qmax=1e3):
        return scale_fns.norm_qrange(self, qmin=qmin, qmax=qmax)

    def norm_max(self, ):
        return scale_fns.norm_max(self)

    
    
    
    

    


    




    








