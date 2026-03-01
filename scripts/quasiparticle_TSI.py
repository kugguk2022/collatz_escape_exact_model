#quasiparticle invariants for type 2 sectors
import numpy as np

def propagation_map(u, v, c, gamma):
    """
    Exactly reversible forward step of the Collatz-Tetration Quasi-Particle.
    """
    # Step 1: Update 'u' using current 'v'
    F_v = gamma * (c ** v) * np.sin(np.pi * v)
    u_next = u + F_v
    
    # Step 2: Update 'v' using the NEW 'u'
    G_u = gamma * (c ** u_next) * np.cos(np.pi * u_next)
    v_next = v + G_u
    
    return u_next, v_next

def retraction_map(u_next, v_next, c, gamma):
    """
    The flawless time-reversed step. Bidirectional arrow.
    """
    # Step 1: Recover 'v' using the known 'u_next'
    G_u = gamma * (c ** u_next) * np.cos(np.pi * u_next)
    v = v_next - G_u
    
    # Step 2: Recover 'u' using the recovered 'v'
    F_v = gamma * (c ** v) * np.sin(np.pi * v)
    u = u_next - F_v
    
    return u, v

def compute_quasi_particle_state(u0, v0, c, gamma):
    """
    Evolves the particle and computes the local Braid Phase (Trace).
    """
    # Evolve forward
    u1, v1 = propagation_map(u0, v0, c, gamma)
    
    # To find the Trace of the Jacobian, we need the derivatives of the gates
    # dF/dv of gamma * c^v * sin(pi*v)
    dF_dv = gamma * (c ** v0) * (np.log(c) * np.sin(np.pi * v0) + np.pi * np.cos(np.pi * v0))
    
    # dG/du of gamma * c^u_next * cos(pi*u_next)
    dG_du = gamma * (c ** u1) * (np.log(c) * np.cos(np.pi * u1) - np.pi * np.sin(np.pi * u1))
    
    # The Trace of this specific symplectic composition is exactly:
    trace_J = 2 + (dF_dv * dG_du)
    
    # Retract back to prove massless bidirectional time arrow
    u_rec, v_rec = retraction_map(u1, v1, c, gamma)
    distance_error = np.abs(u_rec - u0) + np.abs(v_rec - v0)
    
    return trace_J, distance_error

# --- Sweep Example ---
if __name__ == "__main__":
    c_base = complex(1.11, 0)
    
    # A single quasi-particle at a starting coordinate
    u_start, v_start = complex(0.5, 0.1), complex(0.5, 0.1)
    
    print("Gamma  | Distance Error | |Trace| | Regime")
    print("-" * 50)
    
    for gamma_test in [0.01, 0.1, 0.3, 0.5, 0.8]:
        trace, error = compute_quasi_particle_state(u_start, v_start, c_base, gamma_test)
        abs_trace = np.abs(trace)
        
        regime = "GLIDING (Elliptic)" if abs_trace < 2 else "SCATTERING (Hyperbolic)"
        print(f"{gamma_test:4.2f}   | {error:.6e}   | {abs_trace:.4f}  | {regime}")
