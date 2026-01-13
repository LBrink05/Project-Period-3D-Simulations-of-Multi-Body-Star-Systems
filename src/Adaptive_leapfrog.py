import numpy as np
from numba import njit
from pathlib import Path

# Symplectic integrators (leapfrog) preserve energy/phase better for long-term dynamics

# CONSTANTS
NUM_BODIES = 3

# setting cwd directory
CWDDIR = Path.cwd()

def Simulate(data_list, precision, duration):
    mass = np.array([data_list[1], data_list[4], data_list[7]], dtype=np.float64)
    start_pos = np.array([
        [data_list[0][0], data_list[0][1], data_list[0][2]],
        [data_list[3][0], data_list[3][1], data_list[3][2]],
        [data_list[6][0], data_list[6][1], data_list[6][2]],
    ], dtype=np.float64)
    start_vel = np.array([
        [data_list[2][0], data_list[2][1], data_list[2][2]],
        [data_list[5][0], data_list[5][1], data_list[5][2]],
        [data_list[8][0], data_list[8][1], data_list[8][2]],
    ], dtype=np.float64)

    tend = float(duration) * 24.0

    out_dt = 1.0

    dt_max = float(precision)

    dt_min = dt_max/200.0

    eta = 0.01

    frames = position_adaptive(tend, out_dt, NUM_BODIES, start_pos, start_vel, mass, eta, dt_min, dt_max)

    # save CSV files
    out_dir = Path(str(CWDDIR)) / "Simulated_Data"
    out_dir.mkdir(parents=True, exist_ok=True)
    for body in range(NUM_BODIES):
        path = out_dir / f"body{body}.csv"
        np.savetxt(path, frames[body], delimiter=",")


# --- Core physics and integrator (Numba) ---

@njit
def acceleration_components(prior_pos, body, MASS, NUM_BODIES):
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

            r2 = rx * rx + ry * ry + rz * rz + 0.001**2
            r3 = r2 ** 1.5

            ax += MASS[other] * rx / r3
            ay += MASS[other] * ry / r3
            az += MASS[other] * rz / r3

    return ax, ay, az


@njit
def compute_accelerations(prior_pos, MASS, NUM_BODIES, acc_out):
    """
    Fills acc_out[b, :] in place.
    """
    for b in range(NUM_BODIES):
        ax, ay, az = acceleration_components(prior_pos, b, MASS, NUM_BODIES)
        acc_out[b, 0] = ax
        acc_out[b, 1] = ay
        acc_out[b, 2] = az


@njit
def adaptive_leapfrog_inplace(MASS, TIMESTEP, NUM_BODIES, prior_pos, prior_vel, acc, acc_new):
    """
    One leapfrog (velocity Verlet) step in place.
    Uses preallocated acc/acc_new to avoid per-step allocations.
    """


    # acceleration at current positions
    compute_accelerations(prior_pos, MASS, NUM_BODIES, acc)

    # half-step velocity
    for b in range(NUM_BODIES):
        prior_vel[b, 0] += 0.5 * TIMESTEP * acc[b, 0]
        prior_vel[b, 1] += 0.5 * TIMESTEP * acc[b, 1]
        prior_vel[b, 2] += 0.5 * TIMESTEP * acc[b, 2]

    # full-step position
    for b in range(NUM_BODIES):
        prior_pos[b, 0] += prior_vel[b, 0] * TIMESTEP
        prior_pos[b, 1] += prior_vel[b, 1] * TIMESTEP
        prior_pos[b, 2] += prior_vel[b, 2] * TIMESTEP

    # acceleration at updated positions
    compute_accelerations(prior_pos, MASS, NUM_BODIES, acc_new)

    # second half-step velocity
    for b in range(NUM_BODIES):
        prior_vel[b, 0] += 0.5 * TIMESTEP * acc_new[b, 0]
        prior_vel[b, 1] += 0.5 * TIMESTEP * acc_new[b, 1]
        prior_vel[b, 2] += 0.5 * TIMESTEP * acc_new[b, 2]


@njit
def position_adaptive(tend, out_dt, NUM_BODIES, START_POS, START_VEL, MASS, eta, dt_min, dt_max):
    """
    Runs the integrator for TOTAL_STEPS and stores only every SAMPLE_EVERY step.
    Returns frames shaped (NUM_BODIES, OUT_STEPS, 6) with [x,y,z,vx,vy,vz].
    """
    max_out = int(tend * out_dt)
    frames = np.zeros((NUM_BODIES, max_out, 6), dtype=np.float64)

    prior_pos = START_POS.copy()
    prior_vel = START_VEL.copy()

    acc = np.zeros((NUM_BODIES, 3), dtype=np.float64)
    acc_new = np.zeros((NUM_BODIES, 3), dtype=np.float64)

    out_i = 0
    # store t=0
    for b in range(NUM_BODIES):
        frames[b, out_i, 0] = prior_pos[b, 0]
        frames[b, out_i, 1] = prior_pos[b, 1]
        frames[b, out_i, 2] = prior_pos[b, 2]
        frames[b, out_i, 3] = prior_vel[b, 0]
        frames[b, out_i, 4] = prior_vel[b, 1]
        frames[b, out_i, 5] = prior_vel[b, 2]

    out_i += 1
    t = 0
    next_out = out_dt

    Mtot = MASS[0]+MASS[1]+MASS[2]

    while t < tend:
        rmin = min_pair_distance(prior_pos)

        dt = eta * ((rmin ** 3)/Mtot) ** 0.5

        if dt < dt_min:
            dt = dt_min
        elif dt > dt_max:
            dt = dt_max

        if t + dt > tend:
            dt = tend - t

        if t + dt > next_out:
            dt = next_out - t


        adaptive_leapfrog_inplace(MASS, dt, NUM_BODIES, prior_pos, prior_vel, acc, acc_new)
        t += dt


        if t >= next_out - 1e-15:
            for b in range(NUM_BODIES):
                frames[b, out_i, 0] = prior_pos[b, 0]
                frames[b, out_i, 1] = prior_pos[b, 1]
                frames[b, out_i, 2] = prior_pos[b, 2]
                frames[b, out_i, 3] = prior_vel[b, 0]
                frames[b, out_i, 4] = prior_vel[b, 1]
                frames[b, out_i, 5] = prior_vel[b, 2]
            out_i += 1
            next_out += out_dt

            if out_i > max_out:
                break

    return frames[:, :out_i, :]

@njit
def min_pair_distance(pos):

    #pair 0&1
    dx = pos[1,0] - pos[0,0]
    dy = pos[1,1] - pos[0,1]
    dz = pos[1,2] - pos[0,2]
    d01= (dx*dx + dy*dy + dz*dz) ** 0.5

    #pair 0&2
    dx = pos[2,0] - pos[0,0]
    dy = pos[2,1] - pos[0,1]
    dz = pos[2,2] - pos[0,2]
    d02= (dx*dx + dy*dy + dz*dz) ** 0.5

    #pair 1&2
    dx = pos[2,0] - pos[1,0]
    dy = pos[2,1] - pos[1,1]
    dz = pos[2,2] - pos[1,2]
    d12= (dx*dx + dy*dy + dz*dz) ** 0.5

    return min(d01, d02, d12)