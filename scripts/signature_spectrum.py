import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def compute_radial_spectrum(volume):
    """
    Computes the 3D FFT and radial Power Spectral Density (PSD) of a 3D volume.
    """
    # 1. Compute 3D Fast Fourier Transform
    # Normalize by the total number of voxels to preserve energy scaling
    fft_vol = np.fft.fftn(volume) / volume.size
    
    # 2. Shift the zero-frequency component to the center of the spectrum
    fft_vol_shifted = np.fft.fftshift(fft_vol)
    
    # 3. Calculate the Power Spectrum (magnitude squared)
    power_spectrum = np.abs(fft_vol_shifted)**2
    
    # 4. Generate the 3D wavenumber grid (k-space)
    nx, ny, nz = volume.shape
    kx = np.fft.fftshift(np.fft.fftfreq(nx))
    ky = np.fft.fftshift(np.fft.fftfreq(ny))
    kz = np.fft.fftshift(np.fft.fftfreq(nz))
    
    Kx, Ky, Kz = np.meshgrid(kx, ky, kz, indexing='ij')
    
    # Calculate the scalar wavenumber magnitude |k| for every voxel in k-space
    k_magnitude = np.sqrt(Kx**2 + Ky**2 + Kz**2)
    
    # 5. Radial Binning (Isotropic Averaging)
    # Flatten arrays for 1D binning operations
    k_mag_flat = k_magnitude.flatten()
    power_flat = power_spectrum.flatten()
    
    # Define the number of bins up to the Nyquist frequency (0.5 for normalized frequencies)
    num_bins = int(np.min(volume.shape) / 2)
    
    # Compute the mean power within each radial bin using binned_statistic
    radial_mean, bin_edges, _ = stats.binned_statistic(
        k_mag_flat, 
        power_flat, 
        statistic='mean', 
        bins=num_bins, 
        range=(0, 0.5)
    )
    
    # Calculate the center of each frequency bin for plotting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    return bin_centers, radial_mean

# --- Example Execution ---
# Assuming 'collatz_core' is your 3D numpy array containing the decohered phases
# (We simulate a random core here for demonstration)
collatz_core = np.random.choice([-np.pi, np.pi], size=(64, 64, 64)) 

k_radial, psd_radial = compute_radial_spectrum(collatz_core)

# 6. Visualization
plt.figure(figsize=(8, 5))
plt.plot(k_radial, psd_radial, color='indigo', linewidth=2)
plt.yscale('log') # Log scale is standard for spectral density
plt.title('Radial Power Spectral Density of 3D Collatz Core')
plt.xlabel('Spatial Frequency $|\mathbf{k}|$')
plt.ylabel('Power $S(|\mathbf{k}|)$')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.show()