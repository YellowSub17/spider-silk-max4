
import numpy as np
import matplotlib.pyplot as plt
import h5py
from .scan_plot import ScanPlot

from .const import DATA_PATH, PILATUS_PX, EIGER_PX
from .utils import parse_poni


class Scan(ScanPlot):

    def __init__(self, scan_id=102470, det='pilatus'):
        

        self.scan_id = scan_id
        self.raw_path =  f"{DATA_PATH}/raw/scan-{self.scan_id}.h5"
        self.det = det.lower()
        
        assert self.det in ['pilatus', 'eiger'], f'det must be "pilatus" or "eiger". {det=}'
        self.det_px = PILATUS_PX if self.det=='pilatus' else EIGER_PX
        
        self.raw_path = f"{DATA_PATH}/raw/scan-{self.scan_id}.h5"
        self.raw_det_path = f"{DATA_PATH}/raw/scan-{self.scan_id}_{self.det}.h5"
        self.raw_pcap_path = f"{DATA_PATH}/raw/scan-{self.scan_id}_pcap.h5"
        self.azint_path = f"{DATA_PATH}/process/azint/scan-{self.scan_id}_{self.det}_integrated.h5"
        


        #### this only works if the scan was made with the dscan cmd
        #with h5py.File(self.raw_path,'r') as f:
            #self.sample_x = f['/entry/instrument/sams4_x/value'][...]
            #self.n_ims = self.sample_x.size
        

        #with h5py.File(self.raw_det_path,'r') as f:
            #self.eiger_count_time = f['/entry/instrument/eiger/count_time/'][...]
            
        
        with h5py.File(self.azint_path,'r') as f:
            self.poni_str = f['/entry/reduction/input/poni'][...].item().decode("utf-8")
            self.det_center = parse_poni(self.poni_str)
            self.qs = f['/entry/data/radial_axis'][:]
            self.dq = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            self.Is = f['/entry/data/I'][...]
            self.Is_err = f['/entry/data/I_errors'][...]

        self.Is_mean = np.mean(self.Is, axis=0)
            
        with h5py.File(self.raw_path,'r') as f:
            self.n_ims = f[f'/entry/instrument/{self.det}/data/'][...].shape[0]
            


    def get_imgs(self):
        with h5py.File(self.raw_path,'r') as f:
            imgs = f[f'/entry/instrument/{self.det}/data/'][...]
        return imgs
 
        
    



        
            

