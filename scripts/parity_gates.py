import numpy as np

def collatz_parity_gate(x_grid, y_grid, max_int=5000):
    """
    Generates a continuous-bounded spatial parity gate using Collatz stopping times.
    """
    # 1. Normalize the continuous spatial grids to [0, 1]
    x_norm = (x_grid - np.min(x_grid)) / (np.max(x_grid) - np.min(x_grid) + 1e-9)
    y_norm = (y_grid - np.min(y_grid)) / (np.max(y_grid) - np.min(y_grid) + 1e-9)
    
    # 2. Create a heterogeneous "seed" landscape and scale to integers
    # We mix x and y to create nonlinear spatial variations
    spatial_mix = x_norm * np.cos(y_norm * np.pi) + y_norm * np.sin(x_norm * np.pi)
    spatial_mix = (spatial_mix - np.min(spatial_mix)) / (np.max(spatial_mix) - np.min(spatial_mix) + 1e-9)
    
    # Scale to safe integers [1, max_int] to avoid tetration overflow
    grid_int = np.clip(np.int64(spatial_mix * max_int) + 1, 1, max_int)
    
    # 3. Vectorized Collatz Stopping Time Calculation
    stopping_times = np.zeros_like(grid_int, dtype=int)
    current_vals = np.copy(grid_int)
    active_mask = current_vals > 1  # Track which pixels haven't reached 1 yet
    
    max_iters = 1000 # Hard cap to prevent infinite loops
    
    for _ in range(max_iters):
        if not np.any(active_mask):
            break
            
        # Extract the currently active values
        active_vals = current_vals[active_mask]
        
        # Apply Collatz rules: Even -> n/2, Odd -> 3n+1
        even_submask = (active_vals % 2 == 0)
        odd_submask = ~even_submask
        
        active_vals[even_submask] = active_vals[even_submask] // 2
        active_vals[odd_submask] = 3 * active_vals[odd_submask] + 1
        
        # Update the main array and stopping times
        current_vals[active_mask] = active_vals
        stopping_times[active_mask] += 1
        
        # Update which pixels are still active (> 1)
        active_mask = current_vals > 1

    # 4. Map back to [0, 1] for the Kuramoto equation H(Delta; x)
    # Using modulo 2 creates strict 0 or 1 parity regions. 
    # If you want it slightly smoother for the ODE solver, use a cosine map:
    # g_x = 0.5 * (1 - np.cos(np.pi * stopping_times))
    
    g_x = (stopping_times % 2).astype(float)
    
    return g_x