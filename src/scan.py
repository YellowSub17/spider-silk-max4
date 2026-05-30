
import numpy as np
import matplotlib.pyplot as plt
import h5py
from .scan_plot import ScanPlot

from .const import DATA_PATH, PILATUS_PX, EIGER_PX
from .utils import parse_poni


class Scan(ScanPlot):

    def __init__(self, scan_id=102470, det='pilatus', load_imgs=False):
        

        self.scan_id = scan_id
        self.det = det
        assert self.det in ['pilatus', 'eiger'], f'det must be "pilatus" or "eiger". {det=}'
        self.det_px = PILATUS_PX if self.det=='pilatus' else EIGER_PX

        
        
        
        self.raw_path = f"{DATA_PATH}/raw/scan-{self.scan_id}.h5"
        self.raw_det_path = f"{DATA_PATH}/raw/scan-{self.scan_id}_{self.det}.h5"
        self.raw_pcap_path = f"{DATA_PATH}/raw/scan-{self.scan_id}_pcap.h5"
        self.azint_path = f"{DATA_PATH}/process/azint/scan-{self.scan_id}_{self.det}_integrated.h5"
        
        
        with h5py.File(self.azint_path,'r') as f:
            self.poni_str = f['/entry/reduction/input/poni'][...].item().decode("utf-8")
            self.det_center = parse_poni(self.poni_str)
            self.qs = f['/entry/data/radial_axis'][:]
            self.dq = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            self.Is = f['/entry/data/I'][...]
            self.Is_err = f['/entry/data/I_errors'][...]

            
        with h5py.File(self.raw_path,'r') as f:
            self.n_ims = f[f'/entry/instrument/{self.det}/data/'][...].shape[0]
            self.ring_current = f[f'/entry/instrument/ring_current/data'][...]

        with h5py.File(self.raw_pcap_path,'r') as f:
            self.i0 = f['/entry/instrument/pandabox/data/i_0'][...]
            self.it = f['/entry/instrument/pandabox/data/i_t'][...]

        self.flags = (self.i0==self.i0)

        self.imgs=None
        self.img_mean=None
        self.img_std=None
        if load_imgs:
            self.load_img_data()
            

        
        

    def load_img_data(self):
        if self.imgs is not None:
            return
        with h5py.File(self.raw_path,'r') as f:
            imgs = f[f'/entry/instrument/{self.det}/data/'][...]
        self.imgs = imgs
        self.img_mean = imgs.mean(axis=0)
        self.img_std = imgs.std(axis=0)



    
        
        
    



        
            

