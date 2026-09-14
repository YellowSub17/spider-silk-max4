
import numpy as np
import matplotlib.pyplot as plt
import h5py
import hdf5plugin
from scipy import stats


from . import const

from . import plot_fns
from . import scale_fns

from . import det

class Scan:

    def __init__(self, scan_id=102470, load_imgs=True, load_n=1):


        self.scan_id = scan_id


        ### Paths for reading data
        self.raw_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}.h5"
        self.raw_saxs_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}_eiger.h5"
        self.raw_waxs_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}_pilatus.h5"
        self.raw_pcap_path = f"{const.DATA_PATH}/raw/scan-{self.scan_id}_pcap.h5"
        self.azint_saxs_path = f"{const.DATA_PATH}/process/azint/scan-{self.scan_id}_eiger_integrated.h5"
        self.azint_waxs_path = f"{const.DATA_PATH}/process/azint/scan-{self.scan_id}_pilatus_integrated.h5"


        ### Number of exposures in this run. Usually 5, but by default 1 is loaded.
        with h5py.File(self.raw_path,'r') as f:
            if load_n=='all':
                self.n_ims = f[f'/entry/instrument/ring_current/data/'][...].shape[0]
            else:
                self.n_ims = load_n
            self.ring_current = f[f'/entry/instrument/ring_current/data'][0:self.n_ims,...]


        ### Load saxs data
        with h5py.File(self.azint_saxs_path,'r') as f:
            self.qs_saxs = f['/entry/data/radial_axis'][:]
            self.dq_saxs = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            self.Is_saxs = f['/entry/data/I'][0:self.n_ims,...]
            self.Is_err_saxs = f['/entry/data/I_errors'][0:self.n_ims,...]

        ### Load waxs data
        with h5py.File(self.azint_waxs_path,'r') as f:
            self.qs_waxs = f['/entry/data/radial_axis'][:]
            self.dq_waxs = f['/entry/data/radial_axis'][1]-f['/entry/data/radial_axis'][0]
            self.Is_waxs = f['/entry/data/I'][0:self.n_ims,...]
            self.Is_err_waxs = f['/entry/data/I_errors'][0:self.n_ims,...]

        ### Load initital and transmitted intensity
        with h5py.File(self.raw_pcap_path,'r') as f:
            self.i0 = f['/entry/instrument/pandabox/data/i_0'][0:self.n_ims,...]
            self.it = f['/entry/instrument/pandabox/data/i_t'][0:self.n_ims,...]



        #self.qmin_saxs = self.qs_saxs[np.nonzero(self.Is_saxs[0])[0][0]]
        #self.qmin_waxs = self.qs_waxs[np.nonzero(self.Is_waxs[0])[0][0]]
        #self.qmax_saxs = self.qs_saxs[np.nonzero(self.Is_saxs[0])[0][-1]]
        #self.qmax_waxs = self.qs_waxs[np.nonzero(self.Is_waxs[0])[0][-1]]


        #### Create single merged intensity plot

        # Overlap region on eiger and pilatus detectors
        qloc_saxs = (self.qs_saxs >=det.DET_OVERLAP[0]) & (self.qs_saxs <= det.DET_OVERLAP[1])
        qloc_waxs = (self.qs_waxs >=det.DET_OVERLAP[0]) & (self.qs_waxs <= det.DET_OVERLAP[1])
        # scale factor to multiple waxs to match the overlap region
        self.saxswaxs_sf = np.median(self.Is_saxs[:, qloc_saxs], axis=1)/np.median(self.Is_waxs[:, qloc_waxs], axis=1)
        self.saxswaxs_sf = self.saxswaxs_sf.reshape(-1,1)
        #the scaled waxs curve
        Is_waxs_scaled = self.Is_waxs*self.saxswaxs_sf

        # create the array of points that stretch over the whole q range.
        q_overlap = self.qs_saxs[(self.qs_saxs >= det.DET_OVERLAP[0]) & (self.qs_saxs <= det.DET_OVERLAP[1])]
        q_combined = np.concatenate([
            self.qs_saxs[self.qs_saxs < det.DET_OVERLAP[0]],
            q_overlap,
            self.qs_waxs[self.qs_waxs > det.DET_OVERLAP[1]]
        ])
        self.qs = np.unique(q_combined) # Ensure sorted and unique



        # Interpolate the saxs profile
        I_saxs_interp = np.array([np.interp(self.qs, self.qs_saxs, Is_saxs_i) for Is_saxs_i in self.Is_saxs])

        # Interpolate the scaled waxs profile
        I_waxs_scaled = self.Is_waxs * self.saxswaxs_sf
        I_waxs_interp = np.array([np.interp(q_combined, self.qs_waxs, I_waxs_scaled_i) for I_waxs_scaled_i in I_waxs_scaled])

        # Initialise the merged intensity array
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

        # Use a sigmoid function to blend the waxs and saxs
        x_norm = (q_overlap_pts - det.DET_OVERLAP[0]) / (det.DET_OVERLAP[1] - det.DET_OVERLAP[0])
        k = 5
        w_waxs = 1 / (1 + np.exp(-k * (x_norm - 0.2)))  # Centered at 20% into the overlap region
        w_saxs = 1.0 - w_waxs
        I_combined[:, mask_overlap] = (w_saxs * I_saxs_interp[:, mask_overlap]) + (w_waxs * I_waxs_interp[:, mask_overlap])

        # Combined intensity
        self.Is = I_combined



        # Streak determinination
        linfit_q_region = np.where( (self.qs>det.LINEAR_SAXS_RANGE[0]) & (self.qs<det.LINEAR_SAXS_RANGE[1]))[0]

        self.linfit_Is = np.log10(self.Is[:, linfit_q_region])
        self.linfit_qs = np.log10(self.qs[linfit_q_region])

        slopes, intercepts = np.polyfit(self.linfit_qs, self.linfit_Is.T, deg=1)

        # Calculate predicted values to compute R^2 for each row efficiently
        # broadcast log_x across all fitted lines
        y_pred = slopes * self.linfit_qs[:, None] + intercepts  # Shape: (N_points, N_rows)
        y_true = self.linfit_Is.T  # Shape: (N_points, N_rows)

        # Residual sum of squares & Total sum of squares
        ss_res = np.sum((y_true - y_pred) ** 2, axis=0)
        ss_tot = np.sum((y_true - np.mean(y_true, axis=0)) ** 2, axis=0)

        self.r_squared = 1 - (ss_res / ss_tot)



        # initialise imgs varible to None, to say it hasn't been loaded
        self.imgs = None
        if load_imgs:
            self.load_imgs()


    def load_imgs(self):
        # If the imgs has been loaded, dont run
        if self.imgs is not None:
            return

        # Otherwise, load the images
        with h5py.File(self.raw_path,'r') as f:
            self.imgs_saxs = f[f'/entry/instrument/eiger/data/'][0:self.n_ims,...]
            self.imgs_waxs = f[f'/entry/instrument/pilatus/data/'][0:self.n_ims,...]





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






    def norm_qrange(self, qmin=0, qmax=1e3):
        return scale_fns.norm_qrange(self, qmin=qmin, qmax=qmax)

    def norm_max(self,):
        return scale_fns.norm_max(self)

    def norm_i0(self,):
        return scale_fns.norm_i0(self)























