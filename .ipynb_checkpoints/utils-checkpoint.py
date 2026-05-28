

import h5py
import numpy as np
import matplotlib.pyplot as plt
import json



PROP_ID = 20251747
DATA_PATH = f'/data/visitors/cosaxs/{PROP_ID}/2026052808'
RAW_DATA_PATH = f'{DATA_PATH}/raw'
PROC_DATA_PATH = f'{DATA_PATH}/process'
PHOTON_ENERGY = 12400.0
PILATUS_PX = 0.00017199999999999998
EIGER_PX = 0.000075
PILATUS_SHAPE = [1475, 1679]
EIGER_SHAPE = [2162, 2068]



class Scan:

    def __init__(self, scan_id=102465, det='pilatus'):
        

        self.scan_id = scan_id
        self.raw_path =  f"{DATA_PATH}/raw/scan-{self.scan_id}.h5"
        self.det = det.lower()
        assert self.det in ['pilatus', 'eiger'], f'det must be "pilatus" or "eiger". {det=}'
        
        self.raw_eiger_path = f"{DATA_PATH}/raw/scan-{self.scan_id}_eiger.h5"
        self.raw_pilatus_path = f"{DATA_PATH}/raw/scan-{self.scan_id}_pilatus.h5"
        self.raw_pcap_path = f"{DATA_PATH}/raw/scan-{self.scan_id}_pcap.h5"
        self.azint_eiger_path = f"{DATA_PATH}/process/azint/scan-{self.scan_id}_eiger_integrated.h5"
        self.azint_pilatus_path = f"{DATA_PATH}/process/azint/scan-{self.scan_id}_pilatus_integrated.h5"

        with h5py.File(self.raw_path,'r') as f:
            self.sample_x = f['/entry/instrument/sams4_x/value'][...]
            self.det_z = f['/entry/instrument/start_positioners/detector_carriage_z_position'][...]
            self.n_ims = self.sample_x.size
        

        with h5py.File(self.raw_eiger_path,'r') as f:
            pass
            #self.eiger_count_time = f['/entry/instrument/eiger/count_time/'][...]
        
        with h5py.File(self.raw_pilatus_path,'r') as f:
            pass

       
        with h5py.File(self.raw_pcap_path,'r') as f:
            pass
    

        
        with h5py.File(self.azint_eiger_path,'r') as f:
            self.eiger_poni_str = f['/entry/reduction/input/poni'][...].item().decode("utf-8")
            self.eiger_center = _parse_poni(self.eiger_poni_str)
            self.eiger_q = f['/entry/data/radial_axis'][:]
            self.eiger_dq = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            self.eiger_i = f['/entry/data/I'][...]
            self.eiger_err = f['/entry/data/I_errors'][...]
            #self.eiger_poni = _parse_poni(f['/entry/reduction/input/poni'][...])

        
        with h5py.File(self.azint_pilatus_path,'r') as f:
            self.pilatus_poni_str = f['/entry/reduction/input/poni'][...].item().decode("utf-8")
            self.pilatus_center = _parse_poni(self.pilatus_poni_str)
            self.pilatus_q = f['/entry/data/radial_axis'][:]
            self.pilatus_dq = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            self.pilatus_i = f['/entry/data/I'][...]
            self.pilatus_err = f['/entry/data/I_errors'][...]
            #self.pilatus_poni = _parse_poni(f['/entry/reduction/input/poni'][...])





    def get_pilatus_imgs(self):
        with h5py.File(self.raw_pilatus_path,'r') as f:
            imgs = f['/entry/instrument/pilatus/data/'][...]
        return imgs
        
    def get_eiger_imgs(self):
        with h5py.File(self.raw_eiger_path,'r') as f:
            imgs = f['/entry/instrument/eiger/data/'][...]
        return imgs
        
        
    def plot_iq(self):
        ## todo: add option to plot either or both the detectors

        cmap = plt.cm.get_cmap('jet', self.n_ims)
        plt.figure()
        plt.title('Pilatus')
        for i in range(self.n_ims):
            plt.plot(self.pilatus_q,self.pilatus_i[i,:], color=cmap(i))
        plt.xlabel('q [1/A]')
        plt.ylabel('Inten.')

    def plot_pilatus_img(self, i=0):
        imgs = self.get_pilatus_imgs()
        extent_centered = [
            0 - self.pilatus_center[1],          # xmin
            PILATUS_SHAPE[1]*PILATUS_PX - self.pilatus_center[1],    # xmax
            PILATUS_SHAPE[0]*PILATUS_PX - self.pilatus_center[0],   # ymin (bottom)
            0 - self.pilatus_center[0]           # ymax (top)
        ]
        plt.figure()
        plt.imshow(imgs[i], extent=extent_centered)
        plt.colorbar()
        plt.scatter(0,0,marker='+', color='red')

    def plot_eiger_img(self, i=0):
        imgs = self.get_eiger_imgs()
        print(self.eiger_center[1],
              self.eiger_center[0],
              EIGER_SHAPE[1],
              EIGER_PX
             )
              
        extent_centered = [
            0 - self.eiger_center[1],          # xmin
            EIGER_SHAPE[1]*EIGER_PX - self.eiger_center[1],    # xmax
            EIGER_SHAPE[0]*EIGER_PX - self.eiger_center[0],   # ymin (bottom)
            0 - self.eiger_center[0]           # ymax (top)
        ]
        plt.figure()
        plt.imshow(imgs[i], extent=extent_centered)
        plt.colorbar()
        plt.scatter(0,0,marker='+', color='red')
        
        
        
        
def _parse_poni(poni_str):
    raw_data = {}
    for line in poni_str.strip().split("\n"):
       
        if line.startswith("#") or not line:
            continue
        key, value = line.split(":", 1)
        raw_data[key.strip()] = value.strip()
        
    return float(raw_data["Poni1"]), float(raw_data["Poni2"])
    
    


        
            




    
            
            






