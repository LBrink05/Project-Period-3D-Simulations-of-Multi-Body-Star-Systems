import csv
import numpy as np
from scipy.linalg import qr, expm

def read_phase_space(NUM_BODIES, path):
    
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

def read_timestep_sizes(path):
    """
    Read timestep sizes from the timestep_sizes.csv file.
    Returns a numpy array of timestep sizes for each frame.
    """
    timestep_path = str(path) + "/timestep_sizes.csv"
    timestep_sizes = []
    
    with open(timestep_path, "r") as data:
        for line in data:
            value = line.strip()
            if value:
                timestep_sizes.append(float(value))
    
    return np.array(timestep_sizes)

def calculate_phase_space_norm(phase_space_data, frame_idx):
    
    vec = []
    for body_data in phase_space_data:
        if frame_idx < len(body_data):
            vec.extend(body_data[frame_idx])
    return np.linalg.norm(vec)

def calculate_trajectory_error(reference_data, simulated_data, frame_idx):
    
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

def calculate_hamiltonian(phase_space_data, masses, G=1.0, softening = 0.001):

    num_bodies = len(phase_space_data)
    
    # Kinetic energy: T = (1/2) * sum m_i * v_i^2
    T = 0.0
    for i in range(num_bodies):
        pos_vel = phase_space_data[i]
        vx, vy, vz = pos_vel[3], pos_vel[4], pos_vel[5]
        v_squared = vx**2 + vy**2 + vz**2
        T += 0.5 * masses[i] * v_squared
    
    # Potential energy: V = -G * sum sum (m_i * m_j / r_ij)
    V = 0.0
    eps2 = softening * softening
    for i in range(num_bodies):
        for j in range(i + 1, num_bodies):
            xi, yi, zi = phase_space_data[i][0], phase_space_data[i][1], phase_space_data[i][2]
            xj, yj, zj = phase_space_data[j][0], phase_space_data[j][1], phase_space_data[j][2]
            
            dx = xj - xi
            dy = yj - yi
            dz = zj - zi
            r = np.sqrt((dx**2 + dy**2 + dz**2))
            
            if r > 0:  # Avoid division by zero
                V -= G * masses[i] * masses[j] / r
    
    return T + V

def calculate_hamiltonian_error(simulated_data, masses, initial_H, frame_idx, G=1.0):
   
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

    # Calculate initial Hamiltonian from simulated data
    initial_state = [body_data[0] for body_data in simulated_data]
    initial_H = calculate_hamiltonian(initial_state, masses)
    
    def error_at_time(t):

        frame_idx = int(t)
        
        # Calculate trajectory error
        traj_error = calculate_trajectory_error(reference_data, simulated_data, frame_idx)
        
        # Calculate Hamiltonian error
        ham_error = calculate_hamiltonian_error(simulated_data, masses, initial_H, frame_idx)
        
        return traj_error, ham_error
    
    return error_at_time

def calculate_max_error(errors, dt=1.0):
    errors_array = np.array(errors)
    # Filter out infinite values
    finite_errors = errors_array[np.isfinite(errors_array)]
    
    if len(finite_errors) == 0:
        return float('inf')
    
    max_error = finite_errors.max()
    
    return max_error

def compute_jacobian(positions, velocities, masses, G=1.0, softening=0.001):
    """
    Jacobian Matrix Structure    
    J = [      0           I/m     ]
        [ dF/dq            0       ]
    
    """
    num_bodies = len(masses)
    n = 3 * num_bodies  # Total position dimensions
    
    # Initialize Jacobian as 6N x 6N matrix
    J = np.zeros((2*n, 2*n))
    
    # Upper right block: dv/dp = I/m
    for i in range(num_bodies):
        for j in range(3):
            idx = 3*i + j
            J[idx, n + idx] = 1.0 / masses[i]
    
    # Lower left block: dF/dq (force derivatives w.r.t. positions)
    for i in range(num_bodies):
        for j in range(num_bodies):
            if i == j:
                # Diagonal terms: sum of interactions with all other bodies
                for k in range(num_bodies):
                    if k != i:
                        r_vec = positions[k] - positions[i]
                        r = np.linalg.norm(r_vec)
                        r_soft = np.sqrt(r**2 + softening**2)
                        
                        # d(F_i)/d(q_i) contribution from body k
                        for alpha in range(3):
                            for beta in range(3):
                                idx_row = n + 3*i + alpha
                                idx_col = 3*i + beta
                                
                                if alpha == beta:
                                    J[idx_row, idx_col] -= G * masses[k] * (
                                        1.0 / r_soft**3 - 3.0 * r_vec[alpha]**2 / r_soft**5
                                    )
                                else:
                                    J[idx_row, idx_col] -= G * masses[k] * (
                                        -3.0 * r_vec[alpha] * r_vec[beta] / r_soft**5
                                    )
            else:
                # Off-diagonal terms: direct interaction between i and j
                r_vec = positions[j] - positions[i]
                r = np.linalg.norm(r_vec)
                r_soft = np.sqrt(r**2 + softening**2)
                
                for alpha in range(3):
                    for beta in range(3):
                        idx_row = n + 3*i + alpha
                        idx_col = 3*j + beta
                        
                        if alpha == beta:
                            J[idx_row, idx_col] = G * masses[j] * (
                                1.0 / r_soft**3 - 3.0 * r_vec[alpha]**2 / r_soft**5
                            )
                        else:
                            J[idx_row, idx_col] = G * masses[j] * (
                                -3.0 * r_vec[alpha] * r_vec[beta] / r_soft**5
                            )
    
    return J


def calculate_lyapunov_exponents(simulated_data, masses, timestep_sizes, renorm_interval=10, G=1.0, skip_frames=1, verbose=True, use_frame_time=True):
    """
    Calculate Lyapunov exponents using timestep sizes from CSV.
    
    Parameters:
    -----------
    simulated_data : list
        Phase space data for each body
    masses : list
        Mass of each body
    timestep_sizes : np.array
        Array of timestep sizes for each frame
    renorm_interval : int
        Number of processed frames between QR renormalization
    G : float
        Gravitational constant
    skip_frames : int
        Process every nth frame
    verbose : bool
        Print diagnostic information
    use_frame_time : bool
        If True, use frame-based time (frame_idx) for consistent time axis with other plots.
        If False, use cumulative timestep time (actual physical time from integrator).
        Frame-based time is recommended for adaptive timestep integrators to ensure
        the time axis matches the error plots.
    
    Returns:
    --------
    lyapunov_spectrum_sorted : np.array
        Sorted Lyapunov exponents (descending)
    lyapunov_time : float
        Characteristic Lyapunov time (1/lambda_max)
    exponents_over_time : list
        Evolution of exponents at each renormalization step
    time_points : list
        Time points corresponding to exponents_over_time
    sorted_indices : np.array
        Indices mapping sorted exponents to original phase space dimensions
    """
    
    num_bodies = len(simulated_data)
    num_frames = len(simulated_data[0])
    phase_dim = 6 * num_bodies
    
    # Calculate total times for diagnostics
    num_timesteps = min(len(timestep_sizes), num_frames)
    total_timestep_time = np.sum(timestep_sizes[:num_timesteps]) if num_timesteps > 0 else 0
    total_frame_time = num_frames  # Frame-based time (will be converted to seconds by /24 in GUI)
    
    if verbose:
        print(f"Lyapunov calculation: {num_frames} frames, skip={skip_frames}")
        print(f"  Timestep array length: {len(timestep_sizes)}")
        print(f"  Total time from timesteps: {total_timestep_time:.6f} sim units")
        print(f"  Total time from frames: {total_frame_time} frames ({total_frame_time/24:.2f} display seconds)")
        print(f"  Time basis: {'frame-based' if use_frame_time else 'timestep-based'}")
        
        if len(timestep_sizes) > 0:
            print(f"  Timestep range: [{timestep_sizes.min():.2e}, {timestep_sizes.max():.2e}]")
            print(f"  Mean timestep: {timestep_sizes.mean():.2e}")
    
    # Initialize perturbation matrix as identity
    V = np.eye(phase_dim)
    
    # Running sum of logarithms for each exponent
    S = np.zeros(phase_dim)
    
    # Store exponents over time for visualization
    exponents_over_time = []
    time_points = []  # Time points for plotting
    
    renorm_count = 0
    cumulative_timestep_time = 0.0  # Actual physical time from timesteps
    processed = 0
    last_frame = 0
    
    for frame_idx in range(skip_frames, num_frames, skip_frames):
        # Accumulate actual timestep time for all frames since last processed frame
        for i in range(last_frame + 1, min(frame_idx + 1, len(timestep_sizes))):
            cumulative_timestep_time += timestep_sizes[i]
        last_frame = frame_idx
        
        # Determine dt for Jacobian evolution based on time basis choice
        if use_frame_time:
            # Use frame-based time: each frame represents 1 time unit
            # This gives consistent time axis with error plots
            current_time = frame_idx
            if len(time_points) > 0:
                dt = skip_frames  # Time between processed frames
            else:
                dt = frame_idx
        else:
            # Use actual timestep time (physical time from integrator)
            current_time = cumulative_timestep_time
            if len(time_points) > 0:
                dt = current_time - time_points[-1]
            else:
                dt = current_time
        
        # Ensure dt is positive and reasonable
        if dt <= 0:
            dt = skip_frames if use_frame_time else 1e-6
        
        # Extract current state
        positions = np.array([simulated_data[b][frame_idx][:3] for b in range(num_bodies)])
        velocities = np.array([simulated_data[b][frame_idx][3:6] for b in range(num_bodies)])
        
        # Compute Jacobian at current state
        J = compute_jacobian(positions, velocities, masses, G)
        
        # Evolve perturbation vectors
        # Use approximation for small dt, full expm for larger dt
        if dt < 0.1:
            Jdt = J * dt
            V = (np.eye(phase_dim) + Jdt + 0.5 * Jdt @ Jdt) @ V
        else:
            V = expm(J * dt) @ V
        
        processed += 1
        
        # Perform QR decomposition at renormalization intervals
        if (processed % renorm_interval) == 0:
            Q, R = qr(V)
            
            for i in range(phase_dim):
                if abs(R[i, i]) > 0:
                    S[i] += np.log(abs(R[i, i]))
            
            renorm_count += 1
            V = Q
            
            if current_time > 0:
                current_exponents = S / current_time
                exponents_over_time.append(current_exponents.copy())
                time_points.append(current_time)
        
        if verbose and processed % 500 == 0:
            print(f"  Processed {processed}/{num_frames//skip_frames} frames ({100*frame_idx/num_frames:.1f}%)")
    
    if verbose:
        print(f"  Done! Processed {processed} frames, {renorm_count} renormalizations")
        if use_frame_time:
            print(f"  Final time: {time_points[-1] if time_points else 0} frames ({time_points[-1]/24 if time_points else 0:.2f} display seconds)")
        else:
            print(f"  Final time: {time_points[-1] if time_points else 0:.6f} sim units")
    
    # Calculate final Lyapunov spectrum
    total_time = time_points[-1] if time_points else 1.0
    if total_time > 0:
        lyapunov_spectrum = S / total_time
    else:
        lyapunov_spectrum = np.zeros(phase_dim)
    
    # Sort by magnitude (descending)
    sorted_indices = np.argsort(lyapunov_spectrum)[::-1]
    lyapunov_spectrum_sorted = lyapunov_spectrum[sorted_indices]
    
    # Calculate Lyapunov time
    lambda_max = lyapunov_spectrum_sorted[0]
    if lambda_max > 0:
        lyapunov_time = 1.0 / lambda_max
    else:
        lyapunov_time = float('inf')
    
    return lyapunov_spectrum_sorted, lyapunov_time, exponents_over_time, time_points, sorted_indices