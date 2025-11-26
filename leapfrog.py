import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import csv

#symplectic integrators (leapfrog, SABA) preserve energy/phase better for long-term dynamics, while high-order adaptive Runge-Kutta (DOP853) is good for short to medium term with accuracy control

#symplectic integrators
#adaptive timestepping

#Standard ODE solvers (like RK45, DOP853, etc.) integrate accurately in the short term, but they don’t respect the geometry of Hamiltonian mechanics — so they slowly drift in conserved quantities like total energy or angular momentum over long integrations.
#Heuns method (bad)

# CONSTANTS
NUM_BODIES = 3
GRAV = 1 #6.6743015×10E−11 m³/kg*s² Gravitational constant

TIMELINE = np.linspace(0,1000,1000) #length of timeline
MASS = np.array([1,1,1], dtype = float)
TIMESTEP = 0.01 #timestep size (adjust for error) #use variable resolution
TIMESTEP_NUM = int(TIMELINE.size / TIMESTEP) #must be int
FRAMERATIO = int(1 / TIMESTEP) #ratio of frames to timesteps

#function to be called in UI
def Simulate(list):
    MASS = np.array([list[1],list[4],list[7]], dtype = float)
    START_POS = np.array([[list[0][0],list[0][1],list[0][2]],[list[3][0],list[3][1],list[3][2]],[list[6][0],list[6][1],list[6][2]]])
    START_VEL = np.array([[list[2][0],list[2][1],list[2][2]],[list[5][0],list[5][1],list[5][2]],[list[8][0],list[8][1],list[8][2]]])
    pos = position(TIMESTEP, TIMESTEP_NUM, NUM_BODIES, START_POS, START_VEL, MASS)
    frames = pos[:, ::FRAMERATIO, :]

    for body in range(0, NUM_BODIES):
        path = 'Simulated_Data/body' + str(body) + '.csv'
        np.savetxt(path, frames[body], delimiter=",")

#function to calculate acceleration at any configuration of bodies
def acceleration(prior_pos,  body, MASS, NUM_BODIES):
    current_acceleration = np.zeros(3)  
    for other_body in range(NUM_BODIES):
        if other_body != body: 
            r_vec = prior_pos[other_body] - prior_pos[body] # displacement vector
            eps = 0.01
            r_norm = np.linalg.norm(r_vec) 
            current_acceleration += MASS[other_body] * r_vec / (r_norm**2 + eps**2)**1.5
     
    current_acceleration *= GRAV
    return current_acceleration

def classical_leapfrog(pos, body, MASS, TIMESTEP, NUM_BODIES, time ,START_POS,START_VEL, prior_pos, prior_vel):
#symplectic numerical method to approximate chaotic system

    #calculate the velocity and position
    if time == 0:
        #calculate acceleration for one body at one instance
        start_acceleration = np.zeros((NUM_BODIES, 3), dtype=float)
        for body in range(0,NUM_BODIES):
            start_acceleration[body, :]  = acceleration(prior_pos, body, MASS, NUM_BODIES)
        for body in range(0,NUM_BODIES):
            #velocity initialization
            prior_vel[body, :] = START_VEL[body, :] + (1/2) * start_acceleration[body, :] * TIMESTEP
            #update position based on velocity
            prior_pos[body, :] = START_POS[body, :] + prior_vel[body, :] * TIMESTEP

    else:
        #calculate acceleration for one body at one instance
        first_acceleration = np.zeros((NUM_BODIES, 3), dtype=float)
        second_acceleration = np.zeros((NUM_BODIES, 3), dtype=float)

        #update acceleration
        for body in range(0,NUM_BODIES):
            first_acceleration[body, :]  = acceleration(prior_pos, body, MASS, NUM_BODIES)

        #update velocity
        for body in range(0,NUM_BODIES):
            #velocity first half-step
            prior_vel[body, :] += TIMESTEP * 0.5 * first_acceleration[body, :]

        #update position
        for body in range(0,NUM_BODIES):
            #position update based on velocity
            prior_pos[body, :] += prior_vel[body, :] * TIMESTEP

        #update acceletation
        for body in range(0,NUM_BODIES):
            #second acceleration based on new position
            second_acceleration[body, :] = acceleration(prior_pos, body, MASS, NUM_BODIES)
        
        #update velocity
        for body in range(0,NUM_BODIES):
            #velocity second half-step
            prior_vel[body, :] += TIMESTEP * 0.5 * second_acceleration[body, :]

    #returns position of 1 body at 1 moment and prior velocity
    return prior_vel, prior_pos
    

def position(TIMESTEP, TIMESTEP_NUM, NUM_BODIES, START_POS, START_VEL, MASS):
    #function to calculate the position of bodies over time
    #pos is vector pos[body][x,y,z]

    pos = np.empty((NUM_BODIES, TIMESTEP_NUM, 3), dtype=float)

    for body in range(0,NUM_BODIES):
        pos[body][0] = START_POS[body]

    prior_vel = START_VEL
    prior_pos = START_POS

    #for every timestep in the timeline
    for time in range(0, TIMESTEP_NUM):
        #for every body at a moment CHANGE TO BE ONE BIG MATRIX!!!! IT LOSES ENERGY SEQUENTIALLY
         #using leapfrog to calculate position for every moment in time
        prior_vel, prior_pos = classical_leapfrog(pos, body, MASS, TIMESTEP, NUM_BODIES, time ,START_POS,START_VEL, prior_pos, prior_vel)
        pos[:,time] = prior_pos
    return pos