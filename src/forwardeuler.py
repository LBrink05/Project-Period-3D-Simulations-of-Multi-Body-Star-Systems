import numpy as np
from numba import njit
from pathlib import Path

# Simplest possible numerical integrator method (Forward Euler)
NUM_BODIES = 3
CWDDIR = Path.cwd()

def Simulate(data_list, precision, configuration, duration):
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

    timestep = float(precision)
    total_steps = int(duration * 24 / timestep) + 1  # include t=0

    sample_every = max(1, int(1 / timestep))
    out_steps = (total_steps - 1) // sample_every + 1

    frames, timestep_size_list = position_sampled(timestep, total_steps, sample_every, out_steps, NUM_BODIES, start_pos, start_vel, mass)

    out_dir = Path(str(CWDDIR)) / "Simulated_Data" / configuration / "forwardeuler"
    out_dir.mkdir(parents=True, exist_ok=True)
    for body in range(NUM_BODIES):
        path = out_dir / f"body{body}.csv"
        np.savetxt(path, frames[body], delimiter=",")
    
    timestep_path = out_dir / "timestep_sizes.csv"
    np.savetxt(timestep_path, timestep_size_list, delimiter=",")


# Newtonian Gravity (with softening)

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

            r2 = rx * rx + ry * ry + rz * rz + 0.001 **2
            r3 = r2 ** 1.5

            ax += MASS[other] * rx / r3
            ay += MASS[other] * ry / r3
            az += MASS[other] * rz / r3

    return ax, ay, az


@njit
def compute_accelerations(prior_pos, MASS, NUM_BODIES, acc_out):
    for b in range(NUM_BODIES):
        ax, ay, az = acceleration_components(prior_pos, b, MASS, NUM_BODIES)
        acc_out[b, 0] = ax
        acc_out[b, 1] = ay
        acc_out[b, 2] = az


@njit
def forward_euler_step_inplace(MASS, TIMESTEP, NUM_BODIES, prior_pos, prior_vel, acc):
    """
    One Forward Euler step in place using preallocated acc.
    """
    compute_accelerations(prior_pos, MASS, NUM_BODIES, acc)

    for b in range(NUM_BODIES):
        prior_vel[b, 0] += acc[b, 0] * TIMESTEP
        prior_vel[b, 1] += acc[b, 1] * TIMESTEP
        prior_vel[b, 2] += acc[b, 2] * TIMESTEP

        prior_pos[b, 0] += prior_vel[b, 0] * TIMESTEP
        prior_pos[b, 1] += prior_vel[b, 1] * TIMESTEP
        prior_pos[b, 2] += prior_vel[b, 2] * TIMESTEP


@njit
def position_sampled(TIMESTEP, TOTAL_STEPS, SAMPLE_EVERY, OUT_STEPS, NUM_BODIES, START_POS, START_VEL, MASS):
    """
    Runs Forward Euler for TOTAL_STEPS and stores only every SAMPLE_EVERY step.
    Returns frames shaped (NUM_BODIES, OUT_STEPS, 6) with [x,y,z,vx,vy,vz].
    """
    frames = np.zeros((NUM_BODIES, OUT_STEPS, 6), dtype=np.float64)
    timestep_size_list = np.zeros(OUT_STEPS, dtype=np.float64)

    prior_pos = START_POS.copy()
    prior_vel = START_VEL.copy()

    acc = np.zeros((NUM_BODIES, 3), dtype=np.float64)

    out_i = 0
    # store t=0
    for b in range(NUM_BODIES):
        frames[b, out_i, 0] = prior_pos[b, 0]
        frames[b, out_i, 1] = prior_pos[b, 1]
        frames[b, out_i, 2] = prior_pos[b, 2]
        frames[b, out_i, 3] = prior_vel[b, 0]
        frames[b, out_i, 4] = prior_vel[b, 1]
        frames[b, out_i, 5] = prior_vel[b, 2]

    timestep_size_list[out_i] = TIMESTEP * SAMPLE_EVERY
    out_i += 1

    for t in range(1, TOTAL_STEPS):
        forward_euler_step_inplace(MASS, TIMESTEP, NUM_BODIES, prior_pos, prior_vel, acc)

        if t % SAMPLE_EVERY == 0:
            for b in range(NUM_BODIES):
                frames[b, out_i, 0] = prior_pos[b, 0]
                frames[b, out_i, 1] = prior_pos[b, 1]
                frames[b, out_i, 2] = prior_pos[b, 2]
                frames[b, out_i, 3] = prior_vel[b, 0]
                frames[b, out_i, 4] = prior_vel[b, 1]
                frames[b, out_i, 5] = prior_vel[b, 2]
            timestep_size_list[out_i] = TIMESTEP * SAMPLE_EVERY
            out_i += 1

            if out_i >= OUT_STEPS:
                break

    return frames, timestep_size_list
