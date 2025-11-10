from manim import *
import numpy as np
import os
import latex


#disable 3D shading on all Surface-derived objects
_old_surface_init = Surface.__init__
 
def _flat_surface_init(self, *args, **kwargs):
    kwargs["shade_in_3d"] = False   # force lighting off
    return _old_surface_init(self, *args, **kwargs)

Surface.__init__ = _flat_surface_init


#requires to play:
#csc pos file, scene scale, x,y,z minmax

class main_scene(ThreeDScene):
    def construct(self):

        # Load simulation data
        folder = "Simulated Data"
        NUM_BODIES = len(os.listdir(folder))
        TIMELINE = np.linspace(0, 1000, 1000)
        frames = np.empty((NUM_BODIES, TIMELINE.size, 3), dtype=float)
        for body in range(NUM_BODIES):
            path = f"{folder}/body{body}.csv"
            frames[body] = np.genfromtxt(path, delimiter=',')

        # Scene setup

        self.camera.background_color = BLACK
        self.renderer.camera.light_source.move_to(100*IN)

        #these are scene specific constants and have to be conveyed somehow
        SCENE_SCALE = 2.5
        x_min, x_max = -10, 10
        y_min, y_max = -10, 10
        z_min, z_max = -10, 10
        #

        #axes setup
        axes = ThreeDAxes(
            x_range=[x_min, x_max, 1],
            y_range=[y_min, y_max, 1],
            z_range=[z_min, z_max, 1],
            x_length=(x_max - x_min),
            y_length=(y_max - y_min),
            z_length=(z_max - z_min)
        )
        #adds numbers
        axes.add_coordinates()

        self.add(axes)

        # Start rotating the camera for a dynamic 3D effect
        self.begin_ambient_camera_rotation(rate=0.2)
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES, zoom=(1/SCENE_SCALE))
        self.renderer.camera.light_source.move_to(1000 * IN)

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


def star_func(radius, color):
    #core
    core = Circle(radius=radius, fill_color=color, fill_opacity=1, stroke_width=0)
    #gradient for star effect
    halo = Circle(radius=radius+0.01, fill_color=color, fill_opacity=0.25, stroke_width=0)
    halo1 = Circle(radius=radius+0.05, fill_color=color, fill_opacity=0.20, stroke_width=0)
    halo2 = Circle(radius=radius+0.1, fill_color=color, fill_opacity=0.15, stroke_width=0)

    return VGroup(core, halo, halo1, halo2)

def path_func(star, color, frames):
    #path = VMobject()
    #path.set_points_smoothly([*map(lambda p: np.array(p), frames)])
    #path.set_stroke(color=color, width=2, opacity=0.6)

    path = TracedPath(star.get_center, stroke_color=color, stroke_width=2)

    return path