import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import csv
import os, os.path

NUM_BODIES = len(os.listdir('Simulated Data'))

TIMELINE = np.linspace(0,1000,1000) #length of timeline
TIMESTEP = 0.01 #timestep size (adjust for error) #use variable resolution
TIMESTEP_NUM = int(TIMELINE.size / TIMESTEP) #must be int
FRAMERATIO = int(1 / TIMESTEP) #ratio of frames to timesteps

#slice calculated positions to remove timesteps between frames for easy frame by frame position list
frames = np.empty((NUM_BODIES, TIMELINE.size, 3), dtype=float)
for body in range(0,NUM_BODIES):
    path = 'Simulated Data/body'+ str(body) +'.csv'
    frames[body] = np.genfromtxt(path, delimiter=',')

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
    interval=10,
)

axis_dim = 10
ax.set_xlim(-axis_dim, axis_dim)
ax.set_ylim(-axis_dim, axis_dim)
ax.set_zlim(-axis_dim, axis_dim)

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show() 