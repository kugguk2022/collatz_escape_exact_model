import numpy as np
import matplotlib.pyplot as plt

# --- 1. YOUR FUNCTIONS (from previous step) ---
def C_cos_collatz(z):
    cz = np.cos(np.pi * z)
    return (z / 4.0) * (1.0 + cz) + ((3.0 * z + 1.0) / 2.0) * (1.0 - cz)

def orbit_mean_G_dynamic(c, z0, alpha, gamma_base=0.1, N=100, R=80.0):
    z = np.complex128(z0)
    arg_c = np.angle(c)
    ln_abs_c = np.log(np.abs(c)) if np.abs(c) > 0 else -np.inf
    G_sum = 0.0
    steps = 0

    for i in range(1, N + 1): 
        # Dynamic annealing step: gamma grows with time
        gamma_i = gamma_base * (i / 10.0) # Scaled to /10 so it doesn't instantly explode
        
        Cz = C_cos_collatz(z)
        anti_term = gamma_i * np.tanh(np.abs(z)) * np.sin(np.pi * z)
        Cz_kicked = Cz - anti_term

        G_n = (np.imag(Cz_kicked) * (arg_c + alpha * np.exp(np.real(Cz_kicked))) 
               - np.real(Cz_kicked) * ln_abs_c)
        G_sum += float(np.real(G_n))  
        steps += 1

        phi = alpha * np.exp(np.real(Cz_kicked))
        c_eff = c * np.exp(1j * phi)
        z = np.exp(Cz_kicked * np.log(c_eff))

        if (not np.isfinite(z.real)) or (not np.isfinite(z.imag)) or (np.abs(z) > R):
            return (False, G_sum / steps, steps)

    return (True, G_sum / steps, steps)

# --- 2. GENERATING THE DECOHERENCE MAP ---
if __name__ == "__main__":
    # Setup grid matching your previous image bounds
    res = 200 # Resolution (increase for higher quality, lower for faster render)
    re_c = np.linspace(-2.5, 1.5, res)
    im_c = np.linspace(-2.0, 2.0, res)
    
    survival_grid = np.zeros((res, res))
    energy_grid = np.zeros((res, res)) # To track G_mean of survivors
    
    # Parameters
    alpha_val = 0.1  # Set your alpha
    z0_val = 0.5     # Set your starting z
    base_gamma = 0.5 # The comet strength

    print("Running dynamic thermal stress test...")
    
    for row, im in enumerate(im_c):
        for col, re in enumerate(re_c):
            c = re + 1j * im
            # Run the dynamic orbit
            survived, mean_g, steps = orbit_mean_G_dynamic(
                c, z0=z0_val, alpha=alpha_val, gamma_base=base_gamma, N=80
            )
            
            if survived:
                survival_grid[row, col] = 1 # Yellow
                energy_grid[row, col] = mean_g # Track angular momentum
            else:
                survival_grid[row, col] = 0 # Purple Void

    # --- 3. PLOTTING ---
    plt.figure(figsize=(10, 8))
    
    # We use extent to align the axes correctly
    plt.imshow(survival_grid, origin='lower', extent=[-2.5, 1.5, -2.0, 2.0], cmap='viridis')
    
    plt.title(f"Dynamic Decoherence Survivors (Annealing $\gamma$)\nOnly the deep core orbitals survive")
    plt.xlabel("Re(c)")
    plt.ylabel("Im(c)")
    
    plt.colorbar(label="Survival (1 = Survived, 0 = Decohered)")
    plt.show()