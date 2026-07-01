
import numpy as np
from . import const



class Detector:

    def __init__(self, name, px, shape, center, clen, qrange):
        self.name = name
        self.px = px
        self.shape = shape
        self.center = center
        self.clen = clen
        self.qrange = qrange
        self.xrange = ( -center[1], shape[1]*px - center[1])
        self.yrange = ( -center[0], shape[0]*px - center[0])
        


        x_array = np.linspace(self.xrange[0],self.xrange[1], num=self.shape[1])
        y_array = np.linspace(self.yrange[0],self.yrange[1], num=self.shape[0])
        
        self.img_x, self.img_y = np.meshgrid(x_array, y_array)
        self.img_phi = np.arctan2(self.img_y, self.img_x)
        
        self.img_r = np.sqrt( self.img_x**2 + self.img_y**2)
        self.img_q = ((4*np.pi)/const.PHOTON_LAMBDA)*np.sin(0.5*np.arctan2(self.img_r, self.clen))/1e10
        
        self.det_mask =  (np.load(f'{const.DATA_PATH}/{self.name}_mask.npy')==1)
        
        self.q_mask = (self.img_q >= self.qrange[0]) & ( self.img_q <= self.qrange[1])

        self.mask = self.det_mask & self.q_mask
        
        

PILATUS_PX = 0.000172
PILATUS_SHAPE = (1475, 1679)
PILATUS_CENTER = (0.1763777868067733, 0.1929054147362801)
PILATUS_CLEN = 0.4404227673127206
PILATUS_QRANGE = [1.25e-1, 2.125]

EIGER_PX = 0.000075
EIGER_SHAPE = (2162, 2068)
EIGER_CENTER = (0.0734895634276946, 0.007149631494855128)
EIGER_CLEN = 3.4320383606449973
EIGER_QRANGE = [3.1e-3, 3e-1]


DET_OVERLAP = [1.25e-1, 2.5e-1]
FULL_QRANGE = [3.1e-3,2.125]

        

pilatus_det = Detector('pilatus', PILATUS_PX, PILATUS_SHAPE, PILATUS_CENTER, PILATUS_CLEN, PILATUS_QRANGE)
eiger_det = Detector('eiger', EIGER_PX, EIGER_SHAPE, EIGER_CENTER, EIGER_CLEN, EIGER_QRANGE)
