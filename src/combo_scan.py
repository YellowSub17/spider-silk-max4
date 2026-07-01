



from .scan import Scan
from . import const
from . import plot_fns
from . import scale_fns

import numpy as np
import h5py



class ComboScan:


    def __init__(self, scan_ids, load_imgs=False, load_n=1):
        

        self.scan_ids = scan_ids
        #self.det = det
        
        
        self.scans = []
        for i, scan_id in enumerate(self.scan_ids):
            print(f'Loading Scans for ComboScan:\t{i+1}/{len(self.scan_ids)}', end='\r')
            self.scans.append(Scan(scan_id, load_imgs=False, load_n=load_n))
        print('\nDone.')


        

        #self.raw_paths = list(map(lambda scan: scan.raw_path, self.scans))
        #self.raw_det_paths = list(map(lambda scan: scan.raw_det_path, self.scans))
        #self.raw_pcap_paths = list(map(lambda scan: scan.raw_pcap_path, self.scans))
        #self.azint_paths = list(map(lambda scan: scan.azint_path, self.scans))


        self.qs = self.scans[0].qs[:]
        #self.dq = self.scans[0].dq

        
        self.Is = np.concatenate(list(map(lambda scan: scan.Is, self.scans)), axis=0)
        self.Is_err = np.concatenate(list(map(lambda scan: scan.Is_err, self.scans)))

        self.n_ims = np.sum(list(map(lambda scan: scan.n_ims, self.scans)))
        self.scan_n_ims = self.scans[0].n_ims
        self.ring_current = np.concatenate(list(map(lambda scan: scan.ring_current, self.scans)), axis=0)

        
        self.i0 = np.concatenate(list(map(lambda scan: scan.i0, self.scans)), axis=0)
        self.it = np.concatenate(list(map(lambda scan: scan.it, self.scans)), axis=0)

        self.flags = (self.i0==self.i0)

        self.imgs =None
        if load_imgs:
            self.load_imgs()
        
    def load_imgs(self):
        imgs = []
        for i, scan in enumerate(self.scans):
            print(f'Loading images for ComboScan:\t{i+1}/{len(self.scan_ids)}', end='\r')
            scan.load_imgs()
            imgs.append(scan.imgs)
        print('\nDone')
        self.imgs = np.concatenate(imgs, axis=0)
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
    
    

        

        
        


