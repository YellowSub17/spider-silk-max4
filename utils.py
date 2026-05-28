

import h5py
import numpy as np
import matplotlib.pyplot as plt



PROP_ID = 20251747
DATA_PATH = f'/data/visitors/cosaxs/{PROP_ID}/2026052808'
RAW_DATA_PATH = f'{DATA_PATH}/raw'
PROC_DATA_PATH = f'{DATA_PATH}/process'






class Scan:


    def __init__(self, scan_id=102437):

        self.scan_id = scan_id
  
        self.raw_eiger_path = f"{DATA_PATH}/raw/scan-{self.scan_id}_eiger.h5"
        with h5py.File(self.raw_eiger_path,'r') as f:
            self.eiger_img = f['/entry/instrument/eiger/data/'][0]
            self.eiger_count_time = f['/entry/instrument/eiger/count_time/']
            self.eiger_energy = f['entry/instrument/eiger/photon_energy']

        self.raw_pilatus_path = f"{DATA_PATH}/raw/scan-{self.scan_id}_pilatus.h5"
        with h5py.File(self.raw_pilatus_path,'r') as f:
            self.pilatus_img = f['/entry/instrument/pilatus/data/'][0]
            self.pilatus_energy = f['/entry/instrument/pilatus/energy']
            
        self.raw_pcap_path = f"{DATA_PATH}/raw/scan-{self.scan_id}_pcap.h5"
        with h5py.File(self.raw_pcap_path,'r') as f:
            pass
       

        
        self.azint_eiger_path = f"{DATA_PATH}/process/azint/scan-{self.scan_id}_eiger_integrated.h5"
        with h5py.File(self.azint_eiger_path,'r') as f:
            self.eiger_q = f['/entry/data/radial_axis'][:]
            self.eiger_dq = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            self.eiger_i = f['/entry/data/I'][0]
            self.eiger_err = f['/entry/data/I_errors'][0]

        self.azint_pilatus_path = f"{DATA_PATH}/process/azint/scan-{self.scan_id}_pilatus_integrated.h5"
        with h5py.File(self.azint_pilatus_path,'r') as f:
            self.pilatus_q = f['/entry/data/radial_axis'][:]
            self.pilatus_dq = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            self.pilatus_i = f['/entry/data/I'][0]
            self.pilatus_err = f['/entry/data/I_errors'][0]


    
    def plot_iq(self):
        plt.figure()
        plt.plot(self.eiger_q,self.eiger_i)
        plt.figure()
        plt.plot(self.pilatus_q,self.pilatus_i)
        


        
            




    
            
            






