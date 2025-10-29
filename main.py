import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

#symplectic integrators (leapfrog, SABA) preserve energy/phase better for long-term dynamics, while high-order adaptive Runge-Kutta (DOP853) is good for short to medium term with accuracy control

#symplectic integrators
#adaptive timestepping

#Standard ODE solvers (like RK45, DOP853, etc.) integrate accurately in the short term, but they don’t respect the geometry of Hamiltonian mechanics — so they slowly drift in conserved quantities like total energy or angular momentum over long integrations.
#Heuns method (bad)

# CONSTANTS
NUM_BODIES = 3
MASS = np.array([1,1,1], dtype = float)
GRAV = 1 #6.6743015×10−11 m³/kg*s² Gravitational constant

TIMELINE = np.linspace(0,1000,1000) #length of timeline
TIMESTEP = 0.01 #timestep size (adjust for error)
TIMESTEP_NUM = int(TIMELINE.size / TIMESTEP) #must be int
FRAMERATIO = int(1 / TIMESTEP) #ratio of frames to timesteps


#User Interface
print("Please enter which 3 Body Configuration you want to simulate. Select from 0-2")
configuration_name = input()

if configuration_name == "lagrange-triangle-solution" or configuration_name == "triangle-solution" or configuration_name == "0": #requires dt = 0.01
    # Triangle Solution Lagrange
    r = 10.0  # distance from center of mass to each body
    R = r * np.sqrt(3)
    v_mag = np.sqrt(GRAV * MASS[0] / R)

    # 120-degree separation
    START_POS = np.array([
        [ r, 0, 0],
        [-0.5*r,  np.sqrt(3)/2*r, 0],
        [-0.5*r, -np.sqrt(3)/2*r, 0]
    ], dtype=float)

    START_VEL = np.array([
        [0, v_mag, 0],
        [-v_mag*np.sqrt(3)/2, -v_mag/2, 0],
        [ v_mag*np.sqrt(3)/2, -v_mag/2, 0]
    ], dtype=float)

elif configuration_name == "8-solution" or configuration_name == "eight-solution" or configuration_name == "1": #requires dt = 0.01
    # 8 Figure
    D = 2.57429
    START_POS = np.array([
        [D, 0, 0],
        [-D, 0, 0],
        [0, 0, 0]
    ], dtype=float)

    START_VEL = np.array([
        [0.216343, 0.332029,0],
        [0.216343, 0.332029,0],
        [-0.432686, -0.664058,0]
    ], dtype=float)

elif configuration_name == "yingyang" or configuration_name == "2": #requires dt = 0.0001
    # butterfly Figure
    START_POS = np.array([
        [1,0, 0],
        [-0.5,np.sqrt(3/2), 0],
        [-0.5,-np.sqrt(3/2), 0]
    ], dtype=float)

    START_VEL = np.array([
        [0.0,0.5, 0],
        [-0.433,-0.25, 0],
        [0.433,-0.25, 0]
    ], dtype=float)

else:
    print("Invalid configuration. Exiting program.")
    quit()
    
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
    

def position(TIMESTEP, TIMESTEP_NUM, NUM_BODIES, START_POS, START_VEL):
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

# list of all positions over time
pos = position(TIMESTEP,TIMESTEP_NUM, NUM_BODIES, START_POS, START_VEL)

#slice calculated positions to remove timesteps between frames for easy frame by frame position list
frames = pos[:, ::FRAMERATIO, :]

#plotting the data
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

#animation of plot & points
animated_plots = []
points = []

# initialize plots / points on plots
for body in range(0,len(frames)):
    animated_plots.append(ax.plot([],[],[])[0])
    points.append(ax.plot([],[],[], 'ro', markersize=4)[0])
    
def update_data(frame):

    #go through every E.O.M frame by frame (motion has already been calculated)
    for body in range(len(frames)):

        animated_plots[body].set_data(frames[body, :frame, 0], frames[body, :frame, 1])
        animated_plots[body].set_3d_properties(frames[body, :frame, 2])

        points[body].set_data([frames[body, frame, 0]],[frames[body, frame, 1]])
        points[body].set_3d_properties([frames[body, frame, 2]])

    return animated_plots, points

#actual animation function
animation = FuncAnimation(
    fig=fig,
    func=update_data,
    frames=TIMELINE.size,
    interval=1,
)

axis_dim = 10
ax.set_xlim(-axis_dim, axis_dim)
ax.set_ylim(-axis_dim, axis_dim)
ax.set_zlim(-axis_dim, axis_dim)

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show() 