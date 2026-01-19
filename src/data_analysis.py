import csv
import numpy as np
from scipy.linalg import qr, expm

# Softening constant - MUST match the simulation (see acceleration_components)
SOFTENING = 0.001

# Output interval from simulation - frames are saved at fixed intervals of out_dt
OUTPUT_DT = 1.0  # simulation time units between output frames


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
    NOTE: These are the LAST adaptive substep sizes, not the time between output frames.
    The actual time between output frames is fixed at OUTPUT_DT = 1.0
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


def calculate_hamiltonian(phase_space_data, masses, G=1.0, softening=SOFTENING):
    """
    Calculate the Hamiltonian (total energy) of the system.
    Uses Plummer softening: phi = -Gm/sqrt(r^2 + eps^2)
    This MUST match the simulation's force law.
    """
    num_bodies = len(phase_space_data)
    
    # Kinetic energy: T = (1/2) * sum m_i * v_i^2
    T = 0.0
    for i in range(num_bodies):
        pos_vel = phase_space_data[i]
        vx, vy, vz = pos_vel[3], pos_vel[4], pos_vel[5]
        v_squared = vx**2 + vy**2 + vz**2
        T += 0.5 * masses[i] * v_squared
    
    # Potential energy with softening
    V = 0.0
    eps2 = softening * softening
    for i in range(num_bodies):
        for j in range(i + 1, num_bodies):
            xi, yi, zi = phase_space_data[i][0], phase_space_data[i][1], phase_space_data[i][2]
            xj, yj, zj = phase_space_data[j][0], phase_space_data[j][1], phase_space_data[j][2]
            
            dx = xj - xi
            dy = yj - yi
            dz = zj - zi
            
            r_soft = np.sqrt(dx**2 + dy**2 + dz**2 + eps2)
            V -= G * masses[i] * masses[j] / r_soft
    
    return T + V


def calculate_hamiltonian_error(simulated_data, masses, initial_H, frame_idx, G=1.0, softening=SOFTENING):
    current_state = []
    for body_data in simulated_data:
        if frame_idx < len(body_data):
            current_state.append(body_data[frame_idx])
        else:
            return float('inf')
    
    H_current = calculate_hamiltonian(current_state, masses, G, softening)
    delta_H = abs(H_current - initial_H)
    
    if abs(initial_H) < 1e-10:
        return float('inf')
    
    return (delta_H / abs(initial_H)) * 100.0


def error_function(reference_data, simulated_data, masses):
    initial_state = [body_data[0] for body_data in simulated_data]
    initial_H = calculate_hamiltonian(initial_state, masses)
    
    def error_at_time(t):
        frame_idx = int(t)
        traj_error = calculate_trajectory_error(reference_data, simulated_data, frame_idx)
        ham_error = calculate_hamiltonian_error(simulated_data, masses, initial_H, frame_idx)
        return traj_error, ham_error
    
    return error_at_time


def calculate_max_error(errors, dt=1.0):
    errors_array = np.array(errors)
    finite_errors = errors_array[np.isfinite(errors_array)]
    
    if len(finite_errors) == 0:
        return float('inf')
    
    return finite_errors.max()


def compute_jacobian(positions, velocities, masses, G=1.0, softening=SOFTENING):
    """
    Compute the Jacobian of the Hamiltonian equations of motion.
    
    State vector: [q1, q2, ..., qN, v1, v2, ..., vN]
    where qi = (xi, yi, zi) and vi = (vxi, vyi, vzi)
    
    Equations of motion:
        dq_i/dt = v_i
        dv_i/dt = sum_{j≠i} G * m_j * (q_j - q_i) / |q_j - q_i|³_soft
    
    Jacobian structure:
        J = [  0      I   ]   (6N x 6N matrix)
            [  A      0   ]
    
    where A_ij = ∂(dv_i/dt)/∂q_j
    
    Properties for Hamiltonian systems:
        - trace(J) = 0  (Liouville's theorem)
        - Sum of Lyapunov exponents = 0
    """
    num_bodies = len(masses)
    n = 3 * num_bodies  # position dimensions
    
    J = np.zeros((2*n, 2*n))
    
    # Upper right block: ∂(dq/dt)/∂v = I
    for i in range(n):
        J[i, n + i] = 1.0
    
    eps2 = softening * softening
    
    # Lower left block: ∂(dv/dt)/∂q = A
    # a_i = sum_{j≠i} G * m_j * r_ij / |r_ij|³_soft  where r_ij = q_j - q_i
    #
    # ∂a_i/∂q_j = G * m_j * [I/r³ - 3*r*rᵀ/r⁵]  for j ≠ i
    # ∂a_i/∂q_i = -sum_{j≠i} G * m_j * [I/r³ - 3*r*rᵀ/r⁵]
    
    for i in range(num_bodies):
        for j in range(num_bodies):
            if i == j:
                continue
            
            # Vector from body i to body j
            r_vec = positions[j] - positions[i]
            r2 = np.dot(r_vec, r_vec)
            r2_soft = r2 + eps2
            r_soft = np.sqrt(r2_soft)
            r_soft3 = r_soft ** 3
            r_soft5 = r_soft ** 5
            
            # 3x3 derivative block
            for alpha in range(3):
                for beta in range(3):
                    row = n + 3*i + alpha      # acceleration component of body i
                    col_j = 3*j + beta         # position component of body j
                    col_i = 3*i + beta         # position component of body i
                    
                    if alpha == beta:
                        val = G * masses[j] * (1.0/r_soft3 - 3.0*r_vec[alpha]*r_vec[beta]/r_soft5)
                    else:
                        val = G * masses[j] * (-3.0*r_vec[alpha]*r_vec[beta]/r_soft5)
                    
                    J[row, col_j] += val   # ∂a_i/∂q_j
                    J[row, col_i] -= val   # ∂a_i/∂q_i
    
    return J


def calculate_lyapunov_exponents(simulated_data, masses, timestep_sizes, renorm_interval=10, G=1.0, 
                                  softening=SOFTENING, skip_frames=1, verbose=True):
    """
    Calculate Lyapunov exponents from simulation data.
    
    IMPORTANT: The simulation outputs frames at fixed intervals of OUTPUT_DT = 1.0
    simulation time units. The timestep_sizes array contains the last adaptive
    substep size (not the inter-frame time), so we use OUTPUT_DT directly.
    
    Parameters:
    -----------
    simulated_data : list of arrays
        Phase space data [x,y,z,vx,vy,vz] for each body
    masses : list
        Mass of each body
    timestep_sizes : array
        Adaptive substep sizes (for reference only)
    renorm_interval : int
        Frames between QR renormalizations
    G : float
        Gravitational constant
    softening : float
        Plummer softening (must match simulation = 0.001)
    skip_frames : int
        Process every nth frame
    verbose : bool
        Print diagnostics
    
    Returns:
    --------
    lyapunov_spectrum_sorted : sorted exponents (descending)
    lyapunov_time : 1/λ_max
    exponents_over_time : evolution at each renormalization
    time_points : corresponding times (in display seconds = frames/24)
    sorted_indices : index mapping for phase space dimensions
    """
    
    num_bodies = len(simulated_data)
    num_frames = len(simulated_data[0])
    phase_dim = 6 * num_bodies
    
    # Physical time between output frames (fixed by simulation)
    dt_frame = OUTPUT_DT  # = 1.0 simulation time units
    
    # Time step for Jacobian evolution when skipping frames
    dt_physical = dt_frame * skip_frames
    
    if verbose:
        print(f"Lyapunov calculation: {num_frames} frames, skip={skip_frames}")
        print(f"  Output interval: {OUTPUT_DT} sim units")
        print(f"  dt per step: {dt_physical} sim units")
        print(f"  Softening: {softening}")
        print(f"  Total simulation time: {num_frames * OUTPUT_DT} sim units")
    
    # Initialize perturbation matrix (orthonormal basis)
    V = np.eye(phase_dim)
    
    # Running sum of log stretching factors
    S = np.zeros(phase_dim)
    
    # Storage for time evolution
    exponents_over_time = []
    time_points = []
    
    renorm_count = 0
    processed = 0
    
    # Diagnostic: track trace of Jacobian
    trace_values = []
    
    for frame_idx in range(skip_frames, num_frames, skip_frames):
        # Physical time at this frame
        t_physical = frame_idx * OUTPUT_DT
        
        # Extract current state
        positions = np.array([simulated_data[b][frame_idx][:3] for b in range(num_bodies)])
        velocities = np.array([simulated_data[b][frame_idx][3:6] for b in range(num_bodies)])
        
        # Compute Jacobian
        J = compute_jacobian(positions, velocities, masses, G, softening)
        
        # Diagnostic: trace should be 0
        trace_J = np.trace(J)
        trace_values.append(trace_J)
        
        # Evolve perturbation vectors: V' = exp(J * dt) @ V
        Jdt = J * dt_physical
        Jdt_norm = np.linalg.norm(Jdt)
        
        if Jdt_norm < 0.5:
            # Taylor expansion for small arguments
            Jdt2 = Jdt @ Jdt
            Jdt3 = Jdt2 @ Jdt
            expJdt = np.eye(phase_dim) + Jdt + 0.5*Jdt2 + Jdt3/6.0
        else:
            expJdt = expm(Jdt)
        
        V = expJdt @ V
        processed += 1
        
        # QR renormalization
        if processed % renorm_interval == 0:
            Q, R = qr(V)
            
            # Accumulate log of diagonal elements
            for i in range(phase_dim):
                r_ii = abs(R[i, i])
                if r_ii > 1e-300:
                    S[i] += np.log(r_ii)
            
            renorm_count += 1
            V = Q
            
            # Record current exponent estimates
            if t_physical > 0:
                current_exponents = S / t_physical
                exponents_over_time.append(current_exponents.copy())
                # Time in display seconds (frames / 24 fps)
                time_points.append(frame_idx / 24.0)
        
        if verbose and processed % 500 == 0:
            print(f"  Frame {frame_idx}/{num_frames} ({100*frame_idx/num_frames:.1f}%), trace(J)={trace_J:.2e}")
    
    # Final results
    total_physical_time = (num_frames - 1) * OUTPUT_DT
    
    if verbose:
        print(f"  Done: {processed} frames, {renorm_count} renormalizations")
        trace_arr = np.array(trace_values)
        print(f"  trace(J) stats: mean={trace_arr.mean():.2e}, max|trace|={np.abs(trace_arr).max():.2e}")
    
    # Final Lyapunov spectrum
    if total_physical_time > 0:
        lyapunov_spectrum = S / total_physical_time
    else:
        lyapunov_spectrum = np.zeros(phase_dim)
    
    # Sort descending
    sorted_indices = np.argsort(lyapunov_spectrum)[::-1]
    lyapunov_spectrum_sorted = lyapunov_spectrum[sorted_indices]
    
    # Lyapunov time
    lambda_max = lyapunov_spectrum_sorted[0]
    lyapunov_time = 1.0 / lambda_max if lambda_max > 0 else float('inf')
    
    if verbose:
        sum_exp = np.sum(lyapunov_spectrum)
        print(f"  Sum of exponents: {sum_exp:.6e} (should be ≈0 for Hamiltonian)")
        print(f"  λ_max = {lambda_max:.6e}, t_L = {lyapunov_time:.6e}")
    
    return lyapunov_spectrum_sorted, lyapunov_time, exponents_over_time, time_points, sorted_indices


def verify_jacobian_properties(positions, velocities, masses, G=1.0, softening=SOFTENING):
    """
    Verify Jacobian properties for a Hamiltonian system.
    
    Returns dict with:
        - trace: should be exactly 0
        - symplectic_error: J@Ω + Ω@J.T should be 0
        - eigenvalues: for diagnostics
    """
    num_bodies = len(masses)
    n = 3 * num_bodies
    
    J = compute_jacobian(positions, velocities, masses, G, softening)
    
    # 1. Trace
    trace = np.trace(J)
    
    # 2. Symplectic condition: J @ Ω + Ω @ J.T = 0
    Omega = np.zeros((2*n, 2*n))
    Omega[:n, n:] = np.eye(n)
    Omega[n:, :n] = -np.eye(n)
    
    symplectic_test = J @ Omega + Omega @ J.T
    symplectic_error = np.linalg.norm(symplectic_test)
    
    # 3. Eigenvalues
    eigenvalues = np.linalg.eigvals(J)
    
    return {
        'trace': trace,
        'symplectic_error': symplectic_error,
        'eigenvalues': eigenvalues,
        'jacobian': J
    }


def verify_liouville_along_trajectory(simulated_data, masses, G=1.0, softening=SOFTENING, skip_frames=10):
    """
    Verify trace(J) = 0 along the trajectory.
    """
    num_bodies = len(simulated_data)
    num_frames = len(simulated_data[0])
    
    traces = []
    frames = []
    
    for frame_idx in range(0, num_frames, skip_frames):
        positions = np.array([simulated_data[b][frame_idx][:3] for b in range(num_bodies)])
        velocities = np.array([simulated_data[b][frame_idx][3:6] for b in range(num_bodies)])
        
        J = compute_jacobian(positions, velocities, masses, G, softening)
        traces.append(np.trace(J))
        frames.append(frame_idx)
    
    return np.array(traces), np.array(frames)