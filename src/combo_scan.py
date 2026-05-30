



from .scan import Scan
from .scan_plot import ScanPlot
from .const import PILATUS_PX, EIGER_PX
import numpy as np
import h5py



class ComboScan(ScanPlot):


    def __init__(self, scan_ids, det='pilatus', load_imgs=False):
        

        self.scan_ids = scan_ids
        self.det = det
        assert self.det in ['pilatus', 'eiger'], f'det must be "pilatus" or "eiger". {det=}'
        self.det_px = PILATUS_PX if self.det=='pilatus' else EIGER_PX
        
        self.scans = []
        for i, scan_id in enumerate(self.scan_ids):
            print(f'Loading:\t{i}/{len(self.scan_ids)}', end='\r')
            self.scans.append(Scan(scan_id, det=self.det))

        

        self.raw_paths = list(map(lambda scan: scan.raw_path, self.scans))
        self.raw_det_paths = list(map(lambda scan: scan.raw_det_path, self.scans))
        self.raw_pcap_paths = list(map(lambda scan: scan.raw_pcap_path, self.scans))
        self.azint_paths = list(map(lambda scan: scan.azint_path, self.scans))


        self.poni_str = self.scans[0].poni_str
        self.det_center = self.scans[0].det_center
        self.qs = self.scans[0].qs[:]
        self.dq = self.scans[0].dq

        
        self.Is = np.concatenate(list(map(lambda scan: scan.Is, self.scans)), axis=0)
        self.Is_err = np.concatenate(list(map(lambda scan: scan.Is_err, self.scans)))

        self.n_ims = np.sum(list(map(lambda scan: scan.n_ims, self.scans)))
        self.ring_current = np.concatenate(list(map(lambda scan: scan.ring_current, self.scans)), axis=0)

        
        self.i0 = np.concatenate(list(map(lambda scan: scan.i0, self.scans)), axis=0)
        self.it = np.concatenate(list(map(lambda scan: scan.it, self.scans)), axis=0)

        self.flags = (self.i0==self.i0)

        self.imgs=None
        self.img_mean=None
        self.img_std=None
        if load_imgs:
            self.load_img_data()
            

        
        

    def load_img_data(self):
        if self.imgs is not None:
            return
        imgs =[]
        for i, raw_path in enumerate(self.raw_paths):
            print(f'Loading:\t{i}/{len(self.raw_paths)}', end='\r')
            with h5py.File(raw_path,'r') as f:
                imgs.append(f[f'/entry/instrument/{self.det}/data/'][...])
        
        self.imgs = np.concatenate(imgs, axis=0)
        self.img_mean = self.imgs.mean(axis=0)
        self.img_std = self.imgs.std(axis=0)
        
        

        
        


