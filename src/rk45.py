import numpy as np
from pathlib import Path
from scipy.integrate import solve_ivp

# simplest possible numerical integrator method
# CONSTANTS
NUM_BODIES = 3
# setting cwd directory
CWDDIR = Path.cwd()

# function to be called in UI
def Simulate(data_list, precision, duration):
    mass = np.array([data_list[1], data_list[4], data_list[7]], dtype=np.float64)
    start_pos = np.array([[data_list[0][0], data_list[0][1], data_list[0][2]],
                          [data_list[3][0], data_list[3][1], data_list[3][2]],
                          [data_list[6][0], data_list[6][1], data_list[6][2]]], dtype=np.float64)
    start_vel = np.array([[data_list[2][0], data_list[2][1], data_list[2][2]],
                          [data_list[5][0], data_list[5][1], data_list[5][2]],
                          [data_list[8][0], data_list[8][1], data_list[8][2]]], dtype=np.float64)
    TIMELINE = np.linspace(0, duration*24, duration*24)  # length of timeline
    timestep = precision
    timestep_num = int(TIMELINE.size / timestep)  # must be int
    frameratio = int(1 / timestep)  # ratio of frames to timesteps
    pos, vel = position(TIMELINE, NUM_BODIES, start_pos, start_vel, mass)
    frames = np.concatenate([pos[:, ::frameratio, :], vel[:, ::frameratio, :]], axis=-1)
    
    for body in range(0, NUM_BODIES):
        path = Path(str(CWDDIR)) / 'Simulated_Data' / f"body{body}.csv"
        np.savetxt(path, frames[body], delimiter=",")

# function to calculate acceleration at any configuration of bodies
# Newtonian Gravity

def acceleration(prior_pos, body, MASS, NUM_BODIES):
    ax = 0.0
    ay = 0.0
    az = 0.0
    px = prior_pos[body, 0]
    py = prior_pos[body, 1]
    pz = prior_pos[body, 2]
    
    for other in range(NUM_BODIES): 
        if other != body:            
            rx = prior_pos[other, 0] - px
            ry = prior_pos[other, 1] - py
            rz = prior_pos[other, 2] - pz
            r2 = rx * rx + ry * ry + rz * rz + 0.01 * 0.01
            r = r2 ** 1.5
            ax += MASS[other] * rx / r
            ay += MASS[other] * ry / r
            az += MASS[other] * rz / r
    
    return np.array([ax, ay, az], dtype=np.float64)  

def ode_rk45 (t, f, MASS, NUM_BODIES):
    pos = f[:3*NUM_BODIES].reshape(NUM_BODIES, 3)
    vel = f[3*NUM_BODIES:].reshape(NUM_BODIES, 3)
    
    acc = np.zeros_like(pos)

    for b in range(NUM_BODIES):
        acc[b] = acceleration(pos,b, MASS,NUM_BODIES)
    dydt = np.zeros_like(f)
    dydt[:3*NUM_BODIES] =vel.flatten()
    dydt[3*NUM_BODIES:] = acc.flatten()
    return dydt

def position(TIMELINE, NUM_BODIES, START_POS, START_VEL, MASS):
    # function to calculate the position of bodies over time
    # pos is vector pos[body][time_index][x,y,z]
    f0 = np.concatenate([START_POS.flatten(),START_VEL.flatten()])
    rk45 = solve_ivp(fun=ode_rk45,t_span = (TIMELINE[0], TIMELINE[-1]), y0=f0, t_eval = TIMELINE, args = (MASS, NUM_BODIES),method='RK45', rtol=1e-9, atol=1e-12 )
    pos = rk45.y[:3*NUM_BODIES].reshape(NUM_BODIES, 3, -1).transpose(0, 2, 1)
    vel = rk45.y[3*NUM_BODIES:].reshape(NUM_BODIES, 3, -1).transpose(0, 2, 1)
    return pos, vel
