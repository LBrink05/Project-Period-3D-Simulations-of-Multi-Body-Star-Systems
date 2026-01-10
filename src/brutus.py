import numpy as np
from pathlib import Path
import rebound

# simplest possible numerical integrator method
# CONSTANTS
num_bodies = 3
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
    TIMELINE = np.linspace(0, duration * 24, duration * 24)  # length of timeline
    timestep = precision
    timestep_num = int(TIMELINE.size / timestep)  # must be int
    frameratio = int(1 / timestep)  # ratio of frames to timesteps

    pos = position(timestep, timestep_num, num_bodies, start_pos, start_vel, mass)
    frames = pos[:, ::frameratio, :]


    for body in range(0, num_bodies):
        path = Path(str(CWDDIR)) / 'Simulated_Data' / f"body{body}.csv"
        np.savetxt(path, frames[body], delimiter=",")



def position(timestep, timestep_num, num_bodies, start_pos, start_vel, mass):
    # function to calculate the position of bodies over time
    # pos is vector pos[body][x,y,z]

    sim = rebound.Simulation()
    sim.add(m=mass[0], x=start_pos[0][0], y=start_pos[0][1], z=start_pos[0][2], vx=start_vel[0][0], vy=start_vel[0][1],
            vz=start_vel[0][2])
    sim.add(m=mass[1], x=start_pos[1][0], y=start_pos[1][1], z=start_pos[1][2], vx=start_vel[1][0], vy=start_vel[1][1],
            vz=start_vel[1][2])
    sim.add(m=mass[2], x=start_pos[2][0], y=start_pos[2][1], z=start_pos[2][2], vx=start_vel[2][0], vy=start_vel[2][1],
            vz=start_vel[2][2])
    sim.dt = timestep
    sim.status()
    pos = np.zeros((num_bodies, timestep_num, 3), dtype=np.float64)
    # initialize arrays
    prior_pos = start_pos.copy()
    prior_vel = start_vel.copy()
    # store initial pos
    for b in range(num_bodies):
        pos[b, 0] = prior_pos[b]
    # main loop
    for t in range(timestep_num):
        sim.integrate(t*timestep)
        for b in range(num_bodies):
            pos[b,t] = (
                sim.particles[b].x,
                sim.particles[b].y,
                sim.particles[b].z
            )

    return pos