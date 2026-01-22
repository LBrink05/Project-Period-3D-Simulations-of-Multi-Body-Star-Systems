import numpy as np
from pathlib import Path
from scipy.integrate import solve_ivp

# RK45 adaptive integrator - good accuracy, automatic step size control

# CONSTANTS
NUM_BODIES = 3

# setting cwd directory
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

    frames, timestep_size_list = position_sampled(
        timestep, total_steps, sample_every, out_steps,
        NUM_BODIES, start_pos, start_vel, mass
    )

    # save CSV files
    out_dir = Path(str(CWDDIR)) / "Simulated_Data" / configuration / "rk45"
    out_dir.mkdir(parents=True, exist_ok=True)

    for body in range(NUM_BODIES):
        path = out_dir / f"body{body}.csv"
        np.savetxt(path, frames[body], delimiter=",")

    timestep_path = out_dir / "timestep_sizes.csv"
    np.savetxt(timestep_path, timestep_size_list, delimiter=",")


# --- Core physics ---

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


def ode_rk45(t, f, MASS, NUM_BODIES):
    """
    ODE function for scipy's solve_ivp.
    f contains [pos_flat, vel_flat] for all bodies.
    """
    pos = f[:3 * NUM_BODIES].reshape(NUM_BODIES, 3)
    vel = f[3 * NUM_BODIES:].reshape(NUM_BODIES, 3)

    acc = np.zeros((NUM_BODIES, 3), dtype=np.float64)
    for b in range(NUM_BODIES):
        ax, ay, az = acceleration_components(pos, b, MASS, NUM_BODIES)
        acc[b, 0] = ax
        acc[b, 1] = ay
        acc[b, 2] = az

    dydt = np.zeros_like(f)
    dydt[:3 * NUM_BODIES] = vel.flatten()
    dydt[3 * NUM_BODIES:] = acc.flatten()

    return dydt


def position_sampled(TIMESTEP, TOTAL_STEPS, SAMPLE_EVERY, OUT_STEPS, NUM_BODIES, START_POS, START_VEL, MASS):
    """
    Runs the RK45 integrator and samples every SAMPLE_EVERY steps.
    Returns frames shaped (NUM_BODIES, OUT_STEPS, 6) with [x,y,z,vx,vy,vz].
    """
    # Build the timeline of all integration points
    t_end = (TOTAL_STEPS - 1) * TIMESTEP
    t_all = np.linspace(0, t_end, TOTAL_STEPS)

    # Indices to sample (every SAMPLE_EVERY steps)
    sample_indices = np.arange(0, TOTAL_STEPS, SAMPLE_EVERY)
    t_eval = t_all[sample_indices]

    # Initial state vector
    f0 = np.concatenate([START_POS.flatten(), START_VEL.flatten()])

    # Run RK45
    result = solve_ivp(
        fun=ode_rk45,
        t_span=(t_all[0], t_all[-1]),
        y0=f0,
        t_eval=t_eval,
        args=(MASS, NUM_BODIES),
        method='RK45',
        rtol=1e-9,
        atol=1e-12
    )

    # Extract positions and velocities
    # result.y shape: (6*NUM_BODIES, len(t_eval))
    pos = result.y[:3 * NUM_BODIES].reshape(NUM_BODIES, 3, -1).transpose(0, 2, 1)
    vel = result.y[3 * NUM_BODIES:].reshape(NUM_BODIES, 3, -1).transpose(0, 2, 1)

    # Build frames array: (NUM_BODIES, OUT_STEPS, 6)
    actual_out_steps = pos.shape[1]
    frames = np.zeros((NUM_BODIES, actual_out_steps, 6), dtype=np.float64)
    frames[:, :, 0:3] = pos
    frames[:, :, 3:6] = vel

    # Build timestep_size_list from actual evaluation times
    timestep_size_list = np.zeros(actual_out_steps, dtype=np.float64)
    timestep_size_list[0] = TIMESTEP * SAMPLE_EVERY  # first frame
    if actual_out_steps > 1:
        timestep_size_list[1:] = np.diff(result.t)

    return frames, timestep_size_list