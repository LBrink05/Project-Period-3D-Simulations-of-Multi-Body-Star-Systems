import csv
import numpy as np
import numdifftools as nd
import lyapynov

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

def calculate_hamiltonian(phase_space_data, masses, G=1.0):

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

'''def calculate_rms_error(errors, dt=1.0):
    errors_array = np.array(errors)
    # Filter out infinite values
    finite_errors = errors_array[np.isfinite(errors_array)]
    
    if len(finite_errors) == 0:
        return float('inf')
    
    T = len(finite_errors) * dt
    sum_squared = np.sum(finite_errors**2) * dt
    
    return np.sqrt(sum_squared / T)'''

def calculate_max_error(errors, dt=1.0):
    errors_array = np.array(errors)
    # Filter out infinite values
    finite_errors = errors_array[np.isfinite(errors_array)]
    
    if len(finite_errors) == 0:
        return float('inf')
    
    max_error = finite_errors.max()
    
    return max_error

def calculate_lyapunov(simulated_data):

    jacobian_matrices = []

    def calculate_jacobian_matrices_simple(simulated_data):

        total_timesteps = len(simulated_data)
        
        def dynamical_system_derivatives(full_state):
            state = full_state.reshape(num_bodies, 6)
            pos = state[:, :3]
            vel = state[:, 3:]
            
            derivatives = np.zeros_like(full_state)
            
            for b in range(num_bodies):
                # Position derivatives
                derivatives[b*6:b*6+3] = vel[b]
                
                # Acceleration/velocity derivatives
                ax, ay, az = 0.0, 0.0, 0.0
                px, py, pz = pos[b]
                
                for other in range(num_bodies):
                    if other != b:
                        rx = pos[other, 0] - px
                        ry = pos[other, 1] - py
                        rz = pos[other, 2] - pz
                        
                        r2 = rx**2 + ry**2 + rz**2 + 0.001**2
                        r3 = r2**1.5
                        
                        ax += mass[other] * rx / r3
                        ay += mass[other] * ry / r3
                        az += mass[other] * rz / r3
                
                derivatives[b*6+3:b*6+6] = [ax, ay, az]
            
            return derivatives
        
        # Calculate for each time step
        for t in range(total_timesteps):
            current_state = simulated_data[:, t, :].flatten()
            jac_func = nd.Jacobian(dynamical_system_derivatives)
            J_cont = jac_func(current_state)
            
            # Discrete-time approximation
            state_dim = len(current_state)
            phi = np.eye(state_dim) + J_cont * timestep
            jacobian_matrices.append(phi)
        
    return jacobian_matrices

    #discrete_system = lyapynov.DiscreteDS(x0, t0, f, jac)
