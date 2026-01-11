import csv
import numpy as np

def read_phase_space(NUM_BODIES, path):
    """Read phase space data for all bodies from CSV files."""
    phase_space_data = []
    
    for b in range(NUM_BODIES):
        path_b = str(path) + f"/body{b}.csv"
        body_data = []
        
        with open(path_b, "r") as data:
            for line in data:
                values = line.strip().split(",")
                if len(values) >= 6:
                    # [x, y, z, vx, vy, vz]
                    body_data.append([
                        float(values[0]), float(values[1]), float(values[2]),
                        float(values[3]), float(values[4]), float(values[5])
                    ])
        
        phase_space_data.append(np.array(body_data))
    
    return phase_space_data

def calculate_phase_space_norm(phase_space_data, frame_idx):
    """Calculate Euclidean norm of entire phase space at given frame."""
    vec = []
    for body_data in phase_space_data:
        if frame_idx < len(body_data):
            vec.extend(body_data[frame_idx])
    return np.linalg.norm(vec)

def calculate_trajectory_error(reference_data, simulated_data, frame_idx):
    """Calculate percentage trajectory error at given frame."""
    ref_vec = []
    sim_vec = []
    
    for ref_body, sim_body in zip(reference_data, simulated_data):
        if frame_idx < len(ref_body) and frame_idx < len(sim_body):
            ref_vec.extend(ref_body[frame_idx])
            sim_vec.extend(sim_body[frame_idx])
    
    ref_vec = np.array(ref_vec)
    sim_vec = np.array(sim_vec)
    
    ref_norm = np.linalg.norm(ref_vec)
    if ref_norm == 0:
        return float('inf')
    
    diff_norm = np.linalg.norm(ref_vec - sim_vec)
    return (diff_norm / ref_norm) * 100.0

def calculate_hamiltonian(phase_space_data, masses, G=1.0):
    """
    Calculate Hamiltonian (total energy) for a given phase space state.
    H = T + V = (1/2)m*v^2 - G*Σ(m_i*m_j/r_ij)
    
    Args:
        phase_space_data: list of arrays, each [x, y, z, vx, vy, vz] for one body
        masses: list of masses for each body
        G: gravitational constant (default 1.0 for dimensionless units)
    """
    num_bodies = len(phase_space_data)
    
    # Kinetic energy: T = (1/2) * Σ m_i * v_i^2
    T = 0.0
    for i in range(num_bodies):
        pos_vel = phase_space_data[i]
        vx, vy, vz = pos_vel[3], pos_vel[4], pos_vel[5]
        v_squared = vx**2 + vy**2 + vz**2
        T += 0.5 * masses[i] * v_squared
    
    # Potential energy: V = -G * Σ Σ (m_i * m_j / r_ij)
    V = 0.0
    for i in range(num_bodies):
        for j in range(i + 1, num_bodies):
            xi, yi, zi = phase_space_data[i][0], phase_space_data[i][1], phase_space_data[i][2]
            xj, yj, zj = phase_space_data[j][0], phase_space_data[j][1], phase_space_data[j][2]
            
            dx = xj - xi
            dy = yj - yi
            dz = zj - zi
            r = np.sqrt(dx**2 + dy**2 + dz**2)
            
            if r > 0:  # Avoid division by zero
                V -= G * masses[i] * masses[j] / r
    
    return T + V

def calculate_hamiltonian_error(simulated_data, masses, initial_H, frame_idx, G=1.0):
    """Calculate percentage Hamiltonian error at given frame."""
    # Extract current state for all bodies at this frame
    current_state = []
    for body_data in simulated_data:
        if frame_idx < len(body_data):
            current_state.append(body_data[frame_idx])
        else:
            return float('inf')
    
    # Calculate current Hamiltonian
    H_current = calculate_hamiltonian(current_state, masses, G)
    
    # Calculate error
    delta_H = abs(H_current - initial_H)
    
    if abs(initial_H) < 1e-10:  # Avoid division by very small numbers
        return float('inf')
    
    return (delta_H / abs(initial_H)) * 100.0

def error_function(reference_data, simulated_data, masses):
    """
    Returns a function that calculates both trajectory and Hamiltonian errors at time t.
    
    Args:
        reference_data: list of numpy arrays for reference phase space
        simulated_data: list of numpy arrays for simulated phase space
        masses: list of masses [m1, m2, m3, ...]
    
    Returns:
        Function that takes time t (frame number) and returns (E_%, E_H%)
    """
    # Calculate initial Hamiltonian from simulated data
    initial_state = [body_data[0] for body_data in simulated_data]
    initial_H = calculate_hamiltonian(initial_state, masses)
    
    def error_at_time(t):
        """Calculate errors at time t (frame index)."""
        frame_idx = int(t)
        
        # Calculate trajectory error
        traj_error = calculate_trajectory_error(reference_data, simulated_data, frame_idx)
        
        # Calculate Hamiltonian error
        ham_error = calculate_hamiltonian_error(simulated_data, masses, initial_H, frame_idx)
        
        return traj_error, ham_error
    
    return error_at_time

def calculate_rms_error(errors, dt=1.0):
    """
    Calculate RMS error from list of errors.
    E_RMS = sqrt((1/T) * Σ E(t)^2 * Δt)
    """
    errors_array = np.array(errors)
    # Filter out infinite values
    finite_errors = errors_array[np.isfinite(errors_array)]
    
    if len(finite_errors) == 0:
        return float('inf')
    
    T = len(finite_errors) * dt
    sum_squared = np.sum(finite_errors**2) * dt
    
    return np.sqrt(sum_squared / T)