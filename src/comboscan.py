



from .scan import Scan
from .scan_plot import ScanPlot
from .const import PILATUS_PX, EIGER_PX
import numpy as np



class ComboScan(ScanPlot):


    def __init__(self, scan_ids, det='pilatus'):

        self.scan_ids = scan_ids
        self.det = det.lower()
        assert self.det in ['pilatus', 'eiger'], f'det must be "pilatus" or "eiger". {det=}'
        self.det_px = PILATUS_PX if self.det=='pilatus' else EIGER_PX

    
        self.scans = []
        for i, scan_id in enumerate(self.scan_ids):
            print(f'{i}/{len(self.scan_ids)}', end='\r')
            self.scans.append(Scan(scan_id, det=self.det))

        self.raw_paths = list(map(lambda scan: scan.raw_path, self.scans))
        self.raw_det_paths = list(map(lambda scan: scan.raw_det_path, self.scans))
        self.azint_paths = list(map(lambda scan: scan.azint_path, self.scans))
        self.Is = np.concatenate(list(map(lambda scan: scan.Is, self.scans)), axis=0)
        self.Is_mean = np.mean(self.Is, axis=0)
        self.qs = self.scans[0].qs[:]
        self.Is_err = np.concatenate(list(map(lambda scan: scan.Is_err, self.scans)))
    

        self.n_ims = np.sum(list(map(lambda scan: scan.n_ims, self.scans)))
        

