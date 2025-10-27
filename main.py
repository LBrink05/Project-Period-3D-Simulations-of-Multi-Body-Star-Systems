import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

#symplectic integrators (leapfrog, SABA) preserve energy/phase better for long-term dynamics, while high-order adaptive Runge-Kutta (DOP853) is good for short to medium term with accuracy control

#symplectic integrators
#adaptive timestepping

#Standard ODE solvers (like RK45, DOP853, etc.) integrate accurately in the short term, but they don’t respect the geometry of Hamiltonian mechanics — so they slowly drift in conserved quantities like total energy or angular momentum over long integrations.
#Heuns method (bad)

# CONSTANTS
#Solar masses

NUM_BODIES = 3

MASS = np.ones(NUM_BODIES)

#6.6743015×10−11 m³/kg*s² Gravitational constant
GRAV = 1

# Starting positions
START_POS = np.array([
    [5, 0, 0],
    [-2.5, 4.33, 0],
    [-2.5, -4.33, 0]
], dtype=float)

# Starting velocities
START_VEL = np.array([
    [0,0,0],
    [0,0,0],
    [0,0,0]
], dtype=float)


# creating position data for each body (3)
TIMELINE = np.linspace(0,10,100)
TIMESTEP = 1 

def acceleration(prior_pos, time, body, MASS, NUM_BODIES):
    current_acceleration = np.zeros(3)  
    for other_body in range(NUM_BODIES):
        if other_body != body: 
            r_vec = prior_pos[other_body] - prior_pos[body] # displacement vector
            eps = 0.1
            r_norm = np.linalg.norm(r_vec)
            current_acceleration += MASS[other_body] * r_vec / (r_norm**2 + eps**2)**1.5
                
    current_acceleration *= GRAV
    return current_acceleration

def classical_leapfrog(pos, body, MASS, TIMESTEP, NUM_BODIES,time,START_POS,START_VEL, prior_pos, prior_vel):
    #calculate acceleration for one body at one instance
    current_acceleration  = acceleration(prior_pos, time, body, MASS, NUM_BODIES)
    #calculate the velocity and position
    if time == 0:
        #velocity initialization
        vel_init = START_VEL[body, :] + (1/2) * current_acceleration * TIMESTEP
        velocity = vel_init + TIMESTEP * current_acceleration
        prior_vel[body, :] = velocity

        #update position based on velocity
        body_pos = START_POS[body, :] + velocity * TIMESTEP
        prior_pos[body, :] = body_pos
    else:
        #velocity calculation
        velocity = prior_vel[body, :] + TIMESTEP * current_acceleration
        prior_vel[body, :] = velocity

        #position calculation
        body_pos = prior_pos[body, :] + velocity * TIMESTEP
        prior_pos[body, :] = body_pos

    #returns position of 1 body at 1 moment and prior velocity
    return body_pos, prior_vel, prior_pos
    
def position(TIMESTEP, TIMELINE, NUM_BODIES, START_POS, START_VEL):
    #pos is vector pos[body][x,y,z]

    pos = np.empty((NUM_BODIES, TIMELINE.size, 3), dtype=float)

    for body in range(0,NUM_BODIES):
        pos[body][0] = START_POS[body]

    prior_vel = START_VEL
    prior_pos = START_POS

    #for every timestep in the timeline
    for time in range(0, TIMELINE.size):
        #for every body at a moment
        for body in range(0,NUM_BODIES):
            #using leapfrog to calculate position for every moment in time
                body_pos, prior_vel, prior_pos = classical_leapfrog(pos, body, MASS, TIMESTEP, NUM_BODIES,time,START_POS,START_VEL, prior_pos, prior_vel)
                pos[body][time] = body_pos
    return pos

pos = position(TIMESTEP, TIMELINE, NUM_BODIES, START_POS, START_VEL)

#plotting the data
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


#animation of plot & points
animated_plots = []
points = []

# initialize plots / points on plots
for body in range(len(pos)):
    animated_plots.append(ax.plot([],[],[])[0])
    points.append(ax.plot([],[],[], 'ro', markersize=4)[0])
    
def update_data(frame):

    #go through every E.O.M frame by frame (motion has already been calculated)
    for body in range(len(pos)):

        animated_plots[body].set_data(pos[body, :frame, 0], pos[body, :frame, 1])
        animated_plots[body].set_3d_properties(pos[body, :frame, 2])

        points[body].set_data([pos[body, frame, 0]],[pos[body, frame, 1]])
        points[body].set_3d_properties([pos[body, frame, 2]])

        #setting the view limits based on max/min values
        min_x = max(pos[0, frame, 0], pos[1, frame, 0], pos[2, frame, 0] )
        max_x = min(pos[0, frame, 0], pos[1, frame, 0], pos[2, frame, 0] )

        min_y = max(pos[0, frame, 1], pos[1, frame, 1], pos[2, frame, 1] )
        max_y = min(pos[0, frame, 1], pos[1, frame, 1], pos[2, frame, 1] )

        min_z = max(pos[0, frame, 2], pos[1, frame, 2], pos[2, frame, 2] )
        max_z = min(pos[0, frame, 2], pos[1, frame, 2], pos[2, frame, 2] )

        '''ax.set_xlim(max_x,min_x)
        ax.set_ylim(max_y,min_y)
        ax.set_zlim(max_z,min_z)'''

    return animated_plots, points

#actual animation function
animation = FuncAnimation(
    fig=fig,
    func=update_data,
    frames=len(TIMELINE),
    interval=100,
)

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()