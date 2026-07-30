



from .scan import Scan
from . import det
from . import const
from . import plot_fns
from . import scale_fns

import numpy as np
import h5py



class ComboScan:


    def __init__(self, scan_ids, load_imgs=False, load_n=1):

        self.scan_ids = scan_ids

        self.scans = []
        for i, scan_id in enumerate(self.scan_ids):
            print(f'Loading Scans for ComboScan:\t{i+1}/{len(self.scan_ids)}', end='\r')
            self.scans.append(Scan(scan_id, load_imgs=False, load_n=load_n))
        print('\nDone.')



        self.qs = self.scans[0].qs[:]
        self.qs_waxs = self.scans[0].qs_waxs[:]
        self.qs_saxs = self.scans[0].qs_saxs[:]


        self.Is = np.concatenate(list(map(lambda scan: scan.Is, self.scans)), axis=0)
        self.Is_saxs = np.concatenate(list(map(lambda scan: scan.Is_saxs, self.scans)), axis=0)
        self.Is_waxs = np.concatenate(list(map(lambda scan: scan.Is_waxs, self.scans)), axis=0)
        
        self.Is_err_waxs = np.concatenate(list(map(lambda scan: scan.Is_err_waxs, self.scans)))
        self.Is_err_saxs = np.concatenate(list(map(lambda scan: scan.Is_err_saxs, self.scans)))

        self.n_ims = np.sum(list(map(lambda scan: scan.n_ims, self.scans)))
        self.scan_n_ims = self.scans[0].n_ims
        self.ring_current = np.concatenate(list(map(lambda scan: scan.ring_current, self.scans)), axis=0)

        self.i0 = np.concatenate(list(map(lambda scan: scan.i0, self.scans)), axis=0)
        self.it = np.concatenate(list(map(lambda scan: scan.it, self.scans)), axis=0)

        self.r_squared = np.concatenate(list(map(lambda scan: scan.r_squared, self.scans)), axis=0)

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


    def plot_iq(self, i=None, **kwargs):
        if i is None:
            return plot_fns.plot_1Ds(self, self.qs, self.Is, xlim=det.FULL_QRANGE,  **kwargs)
        else:
            return plot_fns.plot_1D(self, self.qs, self.Is[i], xlim=det.FULL_QRANGE,  **kwargs)
        

    def plot_saxs(self, i=None, **kwargs):
        if i is None:
            return plot_fns.plot_1Ds(self, self.qs_saxs, self.Is_saxs, xlim=det.EIGER_QRANGE,  **kwargs)
        else:
            return plot_fns.plot_1D(self, self.qs_saxs, self.Is_saxs[i], xlim=det.EIGER_QRANGE,  **kwargs)

            
    def plot_waxs(self, i=None, **kwargs):
        if i is None:
            return plot_fns.plot_1Ds(self, self.qs_waxs, self.Is_waxs, xlim=det.PILATUS_QRANGE,  **kwargs)
        else:
            return plot_fns.plot_1D(self, self.qs_waxs, self.Is_waxs[i], xlim=det.PILATUS_QRANGE,  **kwargs)
        






