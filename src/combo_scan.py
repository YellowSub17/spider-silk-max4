



from .scan import Scan
from . import const
from . import utils

import numpy as np
import h5py



class ComboScan:


    def __init__(self, scan_ids, det='pilatus', load_imgs=False, load_n='all'):
        

        self.scan_ids = scan_ids
        assert det in ['pilatus', 'eiger', 'saxs', 'waxs', 'sax', 'wax', 's', 'w'], f'det must be one of pilatus, eiger, saxs, waxs, sax or wax. {det=}'
        if det in ['saxs', 'sax', 'eiger', 's']:
            self.det = 'eiger'
        else:
            self.det = 'pilatus'
            
        self.det_px = const.PILATUS_PX if self.det=='pilatus' else const.EIGER_PX
        self.det_shape = const.PILATUS_SHAPE if self.det=='pilatus' else const.EIGER_SHAPE
        
        
        self.scans = []
        for i, scan_id in enumerate(self.scan_ids):
            print(f'Loading Scans for ComboScan:\t{i+1}/{len(self.scan_ids)}', end='\r')
            self.scans.append(Scan(scan_id, det=self.det, load_imgs=False, load_n=load_n))
        print('\nDone.')


        

        self.raw_paths = list(map(lambda scan: scan.raw_path, self.scans))
        self.raw_det_paths = list(map(lambda scan: scan.raw_det_path, self.scans))
        self.raw_pcap_paths = list(map(lambda scan: scan.raw_pcap_path, self.scans))
        self.azint_paths = list(map(lambda scan: scan.azint_path, self.scans))


        self.poni_str = self.scans[0].poni_str
        self.det_center = self.scans[0].det_center
        self.det_d = self.scans[0].det_d
        self.det_lam = self.scans[0].det_lam
        self.qs = self.scans[0].qs[:]
        self.dq = self.scans[0].dq

        
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
            



    def plot_iq(self, **kwargs):
        plot_fns.plot_iq(self, **kwargs)

    def plot_img(self, i=0, **kwargs):
        plot_fns._plot_2d(self.imgs[i], det_center=self.det_center,det_px = self.det_px, **kwargs)
        
            
    def load_imgs(self):
        imgs = []
        for scan in self.scans:
            scan.load_imgs()
            imgs.append(scan.imgs)
        self.imgs = np.concatenate(imgs, axis=0)
        plot_fns.make_img_stats(self)
        

        
        


