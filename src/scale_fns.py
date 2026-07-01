import numpy as np


def norm_qrange(self, qmin=0, qmax=1e3):
    
    qloc = (self.qs >=qmin) & (self.qs <= qmax)
    norm_sf = np.sum(self.Is[:, qloc], axis=1)
    self.Is = self.Is / norm_sf.reshape(-1,1)


def norm_max(self):
    norm_sf = np.max(self.Is, axis=1)
    self.Is=self.Is/norm_sf.reshape(-1,1)

 
  

def x(self):
        
    i = 0
    # 1. Define your q boundaries and angular bins
    r_min = 0.01   # Replace with your desired q_min
    r_max = 0.02# Replace with your desired q_max
    phi_min = -np.pi/2
    phi_max = np.pi/2
    n_bins = 180  # Number of angular bins for phi
    
    # 2. Create a mask for the specific q band and combine it with your detector mask
    r_band_mask = (cs.img_r >= r_min) & (cs.img_r <= r_max)
    phi_band_mask = (cs.img_phi >=phi_min) & (cs.img_phi <=phi_max)
    
    
    combined_mask = r_band_mask & (cs.img_mask > 0) & phi_band_mask  # Exclude dead pixels/gaps
    
    # 3. Flatten and extract the valid intensity and phi values
    valid_intensities = cs.imgs[i, combined_mask] 
    valid_phi = cs.img_phi[combined_mask]
    
    # 4. Sum the intensities in each bin using the weights parameter
    hist_sum, bin_edges = np.histogram(valid_phi, bins=n_bins, weights=valid_intensities)
    
    # 5. Count the number of pixels in each bin
    hist_counts, _ = np.histogram(valid_phi, bins=n_bins)
    
    # 6. Compute the mean intensity per bin
    # Using np.divide with 'where' safely handles any bins that contain 0 pixels
    bin_means = np.divide(hist_sum, hist_counts, out=np.zeros_like(hist_sum, dtype=np.float64), where=hist_counts != 0)
    
    # Calculate the centers of the phi bins for plotting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # 7. Plot the result
    plt.figure(figsize=(8, 4))
    plt.plot(bin_centers, bin_means, color='blue', label='Intensity')
    plt.xlabel('Phi (radians or degrees)')
    plt.ylabel('Average Intensity')
    plt.hlines([np.median(bin_means)], xmin=phi_min, xmax=phi_max, color='red', linestyle='dashed', label='Median')
    plt.legend()
    
    fig, axes = plt.subplots(1,2)
    axes[0].imshow(combined_mask)
    axes[1].imshow(np.log10(cs.imgs[i]*cs.img_mask +1))
    
    
        
        
        