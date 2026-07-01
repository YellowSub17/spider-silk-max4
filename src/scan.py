
import numpy as np
import matplotlib.pyplot as plt
import h5py
import hdf5plugin

#from . import scan_plot 

from . import const

from . import plot_fns
from . import scale_fns

#from .det import pilatus_det, eiger_det
from . import det

class Scan:

    def __init__(self, scan_id=102470, load_imgs=True, load_n=1):
        

        self.scan_id = scan_id

        

        self.raw_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}.h5"
        self.raw_saxs_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}_eiger.h5"
        self.raw_waxs_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}_pilatus.h5"
        
        self.raw_pcap_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}_pcap.h5"
        
        self.azint_saxs_path = f"{const.DATA_PATH}/process/azint/scan-{self.scan_id}_eiger_integrated.h5"
        self.azint_waxs_path = f"{const.DATA_PATH}/process/azint/scan-{self.scan_id}_pilatus_integrated.h5"


        with h5py.File(self.raw_path,'r') as f:
            if load_n=='all':
                self.n_ims = f[f'/entry/instrument/ring_current/data/'][...].shape[0]
            else:
                self.n_ims = load_n
            self.ring_current = f[f'/entry/instrument/ring_current/data'][0:self.n_ims,...]

        
        with h5py.File(self.azint_saxs_path,'r') as f:
            
            self.qs_saxs = f['/entry/data/radial_axis'][:]
            self.dq_saxs = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            self.Is_saxs = f['/entry/data/I'][0:self.n_ims,...]
            self.Is_err_saxs = f['/entry/data/I_errors'][0:self.n_ims,...]

    
        with h5py.File(self.azint_waxs_path,'r') as f:
            
            self.qs_waxs = f['/entry/data/radial_axis'][:]
            self.dq_waxs = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            self.Is_waxs = f['/entry/data/I'][0:self.n_ims,...]
            self.Is_err_waxs = f['/entry/data/I_errors'][0:self.n_ims,...]

        
        with h5py.File(self.raw_pcap_path,'r') as f:
            self.i0 = f['/entry/instrument/pandabox/data/i_0'][0:self.n_ims,...]
            self.it = f['/entry/instrument/pandabox/data/i_t'][0:self.n_ims,...]

        

        self.qmin_saxs = self.qs_saxs[np.nonzero(self.Is_saxs[0])[0][0]]
        self.qmin_waxs = self.qs_waxs[np.nonzero(self.Is_waxs[0])[0][0]]
        
        self.qmax_saxs = self.qs_saxs[np.nonzero(self.Is_saxs[0])[0][-1]]
        self.qmax_waxs = self.qs_waxs[np.nonzero(self.Is_waxs[0])[0][-1]]

        #self.qs_waxs = self.qs_waxs[32:]
        #self.Is_waxs = self.Is_waxs[:, 32:]
        #self.Is_err_waxs = self.Is_err_waxs[:, 32:]
        

        qloc_saxs = (self.qs_saxs >=1.25e-1) & (self.qs_saxs <= 2.5e-1)
        qloc_waxs = (self.qs_waxs >=1.25e-1) & (self.qs_waxs <= 2.5e-1)

        self.saxswaxs_sf = np.median(self.Is_saxs[:, qloc_saxs], axis=1)/np.median(self.Is_waxs[:, qloc_waxs], axis=1)
        self.saxswaxs_sf = self.saxswaxs_sf.reshape(-1,1)
        Is_waxs_scaled = self.Is_waxs*self.saxswaxs_sf#.reshape(-1,1)




            # --- 1. Define your boundaries (from your red dashed lines) ---
        #q_start = det.DET_OVERLAP[0] # Where overlap begins
        #q_end = det.DET_OVERLAP[1]    # Where overlap ends
        
        # --- 2. Define a master q-grid for the final stitched curve ---
        # We take SAXS points up to q_start, and WAXS points after q_end.
        # In the middle, we can use either grid (let's use the SAXS density for the overlap).
        q_overlap = self.qs_saxs[(self.qs_saxs >= det.DET_OVERLAP[0]) & (self.qs_saxs <= det.DET_OVERLAP[1])]
        q_combined = np.concatenate([
            self.qs_saxs[self.qs_saxs < det.DET_OVERLAP[0]],
            q_overlap,
            self.qs_waxs[self.qs_waxs > det.DET_OVERLAP[1]]
        ])
        self.qs = np.unique(q_combined) # Ensure sorted and unique

        
        
        # --- 3. Interpolate both curves onto the master q-grid ---
        # This ensures they share the exact same x-coordinates for easy blending.
        I_saxs_interp = np.array([np.interp(self.qs, self.qs_saxs, Is_saxs_i) for Is_saxs_i in self.Is_saxs])
        
        # Use your calculated scale factor on the WAXS curve
        I_waxs_scaled = self.Is_waxs * self.saxswaxs_sf
        I_waxs_interp = np.array([np.interp(q_combined, self.qs_waxs, I_waxs_scaled_i) for I_waxs_scaled_i in I_waxs_scaled])
        
        # --- 4. Blend them together using a linear ramp ---
        I_combined = np.zeros((self.n_ims, self.qs.size))
        
        # Region 1: Pure SAXS
        mask_pure_saxs = q_combined < det.DET_OVERLAP[0]
        I_combined[:, mask_pure_saxs] = I_saxs_interp[:, mask_pure_saxs]
        
        # Region 2: Pure WAXS
        mask_pure_waxs = q_combined > det.DET_OVERLAP[1]
        I_combined[:, mask_pure_waxs] = I_waxs_interp[:, mask_pure_waxs]
        
        # Region 3: Overlap Blend
        mask_overlap = (q_combined >= det.DET_OVERLAP[0]) & (q_combined <= det.DET_OVERLAP[1])
        q_overlap_pts = q_combined[mask_overlap]
        
        # Calculate weights: 1 at q_start (100% SAXS), 0 at q_end (100% WAXS)
        w_saxs = (det.DET_OVERLAP[1] - q_overlap_pts) / (det.DET_OVERLAP[1] - det.DET_OVERLAP[0])
        w_waxs = 1.0 - w_saxs
        
        
        x_norm = (q_overlap_pts - det.DET_OVERLAP[0]) / (det.DET_OVERLAP[1] - det.DET_OVERLAP[0])
        k = 5  
        w_waxs = 1 / (1 + np.exp(-k * (x_norm - 0.2)))  # Centered at 20% into the overlap region
        w_saxs = 1.0 - w_waxs
        
        I_combined[:, mask_overlap] = (w_saxs * I_saxs_interp[:, mask_overlap]) + (w_waxs * I_waxs_interp[:, mask_overlap])

        self.Is = I_combined
  
                
        self.flags = (self.i0==self.i0)

        self.imgs = None
        if load_imgs:
            self.load_imgs()


    def load_imgs(self):
        if self.imgs is not None:
            return
        
        with h5py.File(self.raw_path,'r') as f:
            self.imgs_saxs = f[f'/entry/instrument/eiger/data/'][0:self.n_ims,...]
            self.imgs_waxs = f[f'/entry/instrument/pilatus/data/'][0:self.n_ims,...]
        
        


    
    def plot_img(self,i=0, **kwargs):
        return plot_fns.plot_img(self,i=i,**kwargs)

    def plot_iq(self, i=0, **kwargs):
        return plot_fns.plot_iq(self, i=i,**kwargs)

    def plot_iqs(self, **kwargs):
        return plot_fns.plot_iqs(self, **kwargs)


    
    def norm_qrange(self, qmin=0, qmax=1e3):
        return scale_fns.norm_qrange(self, qmin=qmin, qmax=qmax)

    def norm_max(self, ):
        return scale_fns.norm_max(self)

    
    
    
    

    


    




    








