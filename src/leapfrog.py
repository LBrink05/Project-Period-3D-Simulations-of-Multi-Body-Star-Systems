import numpy as np
from numba import njit
from pathlib import Path

#symplectic integrators (leapfrog, SABA) preserve energy/phase better for long-term dynamics, while high-order adaptive Runge-Kutta (DOP853) is good for short to medium term with accuracy control

#symplectic integrators
#adaptive timestepping

#Standard ODE solvers (like RK45, DOP853, etc.) integrate accurately in the short term, but they don’t respect the geometry of Hamiltonian mechanics — so they slowly drift in conserved quantities like total energy or angular momentum over long integrations.
#Heuns method (bad)

# CONSTANTS
NUM_BODIES = 3
GRAV = 1 #6.6743015×10E−11 m³/kg*s² Gravitational constant

#setting cwd directory
CWDDIR = Path.cwd()

#function to be called in UI
def Simulate(data_list, precision, duration):
    mass = np.array([data_list[1],data_list[4],data_list[7]], dtype=np.float64)
    start_pos = np.array([[data_list[0][0],data_list[0][1],data_list[0][2]],[data_list[3][0],data_list[3][1],data_list[3][2]],[data_list[6][0],data_list[6][1],data_list[6][2]]], dtype=np.float64)
    start_vel = np.array([[data_list[2][0],data_list[2][1],data_list[2][2]],[data_list[5][0],data_list[5][1],data_list[5][2]],[data_list[8][0],data_list[8][1],data_list[8][2]]], dtype=np.float64)

    TIMELINE = np.linspace(0, duration*24, duration*24)  # length of timeline
    timestep = precision
    timestep_num = int(TIMELINE.size / timestep)  # must be int
    frameratio = int(1 / timestep)  # ratio of frames to timesteps

    pos = position(timestep, timestep_num, NUM_BODIES, start_pos, start_vel, mass)
    frames = pos[:, ::frameratio, :]

    for body in range(0, NUM_BODIES):
        path =  Path(str(CWDDIR)) / 'Simulated_Data' / f"body{body}.csv"
        np.savetxt(path, frames[body], delimiter=",")
    
#function to calculate acceleration at any configuration of bodies
@njit
def acceleration(prior_pos,  body, MASS, NUM_BODIES):
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
@njit
def classical_leapfrog(MASS, TIMESTEP, NUM_BODIES, time ,START_POS,START_VEL, prior_pos, prior_vel):
#symplectic numerical method to approximate chaotic system

    #calculate the velocity and position
    if time == 0:
        # First step: half-step velocity
        # acceleration at t0
        acc0 = np.zeros((NUM_BODIES, 3), dtype=np.float64)

        for b in range(NUM_BODIES):
            acc0[b] = acceleration(prior_pos, b, MASS, NUM_BODIES)

        # half-step velocity update
        for b in range(NUM_BODIES):
            prior_vel[b] = START_VEL[b] + 0.5 * acc0[b] * TIMESTEP
            prior_pos[b] = START_POS[b] + prior_vel[b] * TIMESTEP

    else:
        # Standard leapfrog step
        acc1 = np.zeros((NUM_BODIES, 3), dtype=np.float64)
        for b in range(NUM_BODIES):
            acc1[b] = acceleration(prior_pos, b, MASS, NUM_BODIES)

        # half-step velocity
        for b in range(NUM_BODIES):
            prior_vel[b] += 0.5 * TIMESTEP * acc1[b]

        # full-step position
        for b in range(NUM_BODIES):
            prior_pos[b] += prior_vel[b] * TIMESTEP

        # second acceleration
        acc2 = np.zeros((NUM_BODIES, 3), dtype=np.float64)
        for b in range(NUM_BODIES):
            acc2[b] = acceleration(prior_pos, b, MASS, NUM_BODIES)

        # half-step velocity again
        for b in range(NUM_BODIES):
            prior_vel[b] += 0.5 * TIMESTEP * acc2[b]

    return prior_vel, prior_pos
    
@njit
def position(TIMESTEP, TIMESTEP_NUM, NUM_BODIES, START_POS, START_VEL, MASS):
    #function to calculate the position of bodies over time
    #pos is vector pos[body][x,y,z]

    pos = np.zeros((NUM_BODIES, TIMESTEP_NUM, 3), dtype=np.float64)

    # initialize arrays
    prior_pos = START_POS.copy()
    prior_vel = START_VEL.copy()

    # store initial pos
    for b in range(NUM_BODIES):
        pos[b, 0] = prior_pos[b]

    # main loop
    for t in range(TIMESTEP_NUM):
        prior_vel, prior_pos = classical_leapfrog(MASS, TIMESTEP, NUM_BODIES, t, START_POS, START_VEL, prior_pos, prior_vel)
        for b in range(NUM_BODIES):
            pos[b, t] = prior_pos[b]

    return pos
