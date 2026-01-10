import numpy as np
from numba import njit
from pathlib import Path

# symplectic integrators (leapfrog) preserve energy/phase better for long-term dynamics

# CONSTANTS
NUM_BODIES = 3

# setting cwd directory
CWDDIR = Path.cwd()

# function to be called in UI
def Simulate(data_list, precision, duration):
    mass = np.array([data_list[1], data_list[4], data_list[7]], dtype=np.float64)
    start_pos = np.array([
        [data_list[0][0], data_list[0][1], data_list[0][2]],
        [data_list[3][0], data_list[3][1], data_list[3][2]],
        [data_list[6][0], data_list[6][1], data_list[6][2]]
    ], dtype=np.float64)
    start_vel = np.array([
        [data_list[2][0], data_list[2][1], data_list[2][2]],
        [data_list[5][0], data_list[5][1], data_list[5][2]],
        [data_list[8][0], data_list[8][1], data_list[8][2]]
    ], dtype=np.float64)

    # calculate number of steps
    timestep = precision
    TIMESTEP_NUM = int(duration * 24 / timestep) + 1  # include t=0

    # compute positions
    pos, vel = position(timestep, TIMESTEP_NUM, NUM_BODIES, start_pos, start_vel, mass)

    # select frames to save (1 per unit time)
    frameratio = max(1, int(1 / timestep))
    frames = np.concatenate([pos[:, ::frameratio, :], vel[:, ::frameratio, :]], axis=-1)
    # save CSV files
    for body in range(NUM_BODIES):
        path = Path(str(CWDDIR)) / 'Simulated_Data' / f"body{body}.csv"
        path.parent.mkdir(parents=True, exist_ok=True)  # ensure folder exists
        np.savetxt(path, frames[body], delimiter=",")

# function to calculate acceleration at any configuration of bodies
@njit
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

            r2 = rx * rx + ry * ry + rz * rz + 0.001**2  # softening
            r3 = r2 ** 1.5

            ax += MASS[other] * rx / r3
            ay += MASS[other] * ry / r3
            az += MASS[other] * rz / r3

    return np.array([ax, ay, az], dtype=np.float64)

@njit
def classical_leapfrog(MASS, TIMESTEP, NUM_BODIES, time, START_POS, START_VEL, prior_pos, prior_vel):
    # symplectic numerical method to approximate chaotic system

    # compute acceleration at current positions
    acc = np.zeros((NUM_BODIES, 3), dtype=np.float64)
    for b in range(NUM_BODIES):
        acc[b] = acceleration(prior_pos, b, MASS, NUM_BODIES)

    # half-step velocity update
    for b in range(NUM_BODIES):
        prior_vel[b] += 0.5 * TIMESTEP * acc[b]

    # full-step position update
    for b in range(NUM_BODIES):
        prior_pos[b] += prior_vel[b] * TIMESTEP

    # compute new acceleration at updated positions
    acc_new = np.zeros((NUM_BODIES, 3), dtype=np.float64)
    for b in range(NUM_BODIES):
        acc_new[b] = acceleration(prior_pos, b, MASS, NUM_BODIES)

    # complete velocity update with second half-step
    for b in range(NUM_BODIES):
        prior_vel[b] += 0.5 * TIMESTEP * acc_new[b]

    return prior_vel, prior_pos

@njit
def position(TIMESTEP, TIMESTEP_NUM, NUM_BODIES, START_POS, START_VEL, MASS):
    # function to calculate the position of bodies over time
    # pos is vector pos[body][time_index][x,y,z]

    pos = np.zeros((NUM_BODIES, TIMESTEP_NUM, 3), dtype=np.float64)
    vel = np.zeros((NUM_BODIES, TIMESTEP_NUM, 3), dtype=np.float64)
    prior_pos = START_POS.copy()
    prior_vel = START_VEL.copy()

    # store initial positions at t=0
    for b in range(NUM_BODIES):
        pos[b, 0] = prior_pos[b]
        vel[b, 0] = prior_vel[b]

    # main loop over time steps
    for t in range(1, TIMESTEP_NUM):
        prior_vel, prior_pos = classical_leapfrog(MASS, TIMESTEP, NUM_BODIES, t, START_POS, START_VEL, prior_pos, prior_vel)
        for b in range(NUM_BODIES):
            pos[b, t] = prior_pos[b]
            vel[b, t] = prior_vel[b]

    return pos, vel
