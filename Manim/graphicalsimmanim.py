from manim import *
import numpy as np
import os
import latex


#export MGL_BACKEND=pyglet
#export PYOPENGL_PLATFORM=x11
#manim -pqh Manim/graphicalsimmanim.py MainScene

#requires to play:
#csc pos file, scene scale, x,y,z minmax

class MainScene(ThreeDScene):
    def construct(self):
            
        #these are scene specific constants and have to be conveyed somehow
        SCENE_SCALE = 2.5
        x_min, x_max = -10, 10
        y_min, y_max = -10, 10
        z_min, z_max = -10, 10
        #

        # Scene setup
        self.wait(0.1) 
        self.camera.background_color = BLACK

        self.renderer.camera.near = 0.01
        self.renderer.camera.far  = 1000

        # Start rotating the camera for a dynamic 3D effect
        self.begin_ambient_camera_rotation(rate=0.2)
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES, zoom=(1/SCENE_SCALE))
        self.renderer.camera.light_source.move_to(1000 * IN)

        #grid setup
        # Scale and rotate grid to align with XY plane in 3D space
        grid = grid_func(x_min, x_max, y_min, y_max)
        self.add(grid)

        # Load simulation data
        folder = "Simulated Data"
        NUM_BODIES = len(os.listdir(folder))
        TIMELINE = np.linspace(0, 1000, 1000)
        frames = np.empty((NUM_BODIES, TIMELINE.size, 3), dtype=float)
        for body in range(NUM_BODIES):
            path = f"{folder}/body{body}.csv"
            frames[body] = np.genfromtxt(path, delimiter=',')

        #axes setup
        axes = ThreeDAxes(
            x_range=[x_min, x_max, 1],
            y_range=[y_min, y_max, 1],
            z_range=[z_min, z_max, 1],
            x_length=(x_max - x_min -2),
            y_length=(y_max - y_min -2),
            z_length=(z_max - z_min -2)
        )
        #adds numbers
        axes.add_coordinates()
        self.add(axes)

        #setting up objects in scene
        colors = [RED,GREEN,BLUE]

        #initializing list of objects
        stars = []
        paths = []

        for body in range(0,NUM_BODIES):
            color = colors[body]
            star = star_func(radius=0.15, color=color)

            # Add as fixed-orientation 3D object (always faces camera)
            self.add_fixed_orientation_mobjects(star)
            star.move_to(frames[body][0] * SCENE_SCALE)

            stars.append(star)
            self.add(stars[body])

            path = path_func(star, color=color, frames=frames[body])
            self.add(path)

        #animate their positions 
        def update_all(mob,alpha):
            frame_index = int(alpha * (len(TIMELINE) - 1))
            for body, sphere in enumerate(stars):
                sphere.move_to(frames[body][frame_index])

        self.play(UpdateFromAlphaFunc(VGroup(*stars), update_all), run_time=10, rate_func=linear)
        
        # Stop camera rotation
        self.stop_ambient_camera_rotation()

#function that draws star
def star_func(radius, color):
    #core
    core = Circle(radius=radius, fill_color=color, fill_opacity=1, stroke_width=0)
    #gradient for star effect
    halo = Circle(radius=radius+0.01, fill_color=color, fill_opacity=0.25, stroke_width=0)
    halo1 = Circle(radius=radius+0.05, fill_color=color, fill_opacity=0.20, stroke_width=0)
    halo2 = Circle(radius=radius+0.1, fill_color=color, fill_opacity=0.15, stroke_width=0)

    return VGroup(core, halo, halo1, halo2)

#function to draw path taken by stars
def path_func(star, color, frames):
    #path = VMobject()
    #path.set_points_smoothly([*map(lambda p: np.array(p), frames)])
    #path.set_stroke(color=color, width=2, opacity=0.6)

    path = TracedPath(star.get_center, stroke_color=color, stroke_width=2)

    return path

#function that creates grid
def grid_func(x_min, x_max, y_min, y_max):
    # XY grid plane
    # Vertical grid lines
    step = 0.5
    thickness = 0.001
    stroke_opacity = 0.05
    grid = VGroup()
    grid.z_index = -10
     # vertical (constant x) lines
    for x in np.arange(x_min, x_max + step, step):
        start = np.array([x, y_min, surface_z(x, y_min)])
        end   = np.array([x, y_max, surface_z(x, y_max)])
        grid.add(Line3D(start=start, end=end, thickness=thickness, stroke_color=GREY, stroke_opacity=stroke_opacity))

    # horizontal (constant y) lines
    for y in np.arange(y_min, y_max + step, step):
        start = np.array([x_min, y, surface_z(x_min, y)])
        end   = np.array([x_max, y, surface_z(x_max, y)])
        grid.add(Line3D(start=start, end=end, thickness=thickness, stroke_color=GREY, stroke_opacity=stroke_opacity))

    return grid

#Surface function for grid (used later for gravity visualization)
def surface_z(x, y):
    """Define the curvature of the grid (z = f(x, y))"""
    return 0