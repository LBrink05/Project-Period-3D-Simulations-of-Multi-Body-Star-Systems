from manim import *
import numpy as np
import os
from pathlib import Path


# export MGL_BACKEND=pyglet
# export PYOPENGL_PLATFORM=x11
# manim -pqh --fps 24 Manim/graphicalsimmanim.py MainScene


class MainScene(ThreeDScene):
    def construct(self):
        # Set frame rate
        config.frame_rate = 24
        
        # Scene setup
        self.wait(0.1)
        self.camera.background_color = BLACK
        self.renderer.camera.near = 0.01
        self.renderer.camera.far = 1000

        # Load simulation data first to determine bounds
        script_dir = Path(__file__).parent.resolve()
        folder = script_dir / ".." / "src" / "Simulated_Data"

        NUM_BODIES = len(os.listdir(folder)) - 1

        # Load first file to get dimensions
        sample = np.genfromtxt(folder / "body0.csv", delimiter=',')
        num_timesteps = sample.shape[0]

        # Load all position data
        frames = np.empty((NUM_BODIES, num_timesteps, 3), dtype=float)
        for body in range(NUM_BODIES):
            path = folder / f"body{body}.csv"
            data = np.genfromtxt(path, delimiter=',')
            frames[body] = data[:, :3]

        # Calculate actual bounds from simulation data with padding
        all_positions = frames.reshape(-1, 3)
        x_min_data, y_min_data, z_min_data = all_positions.min(axis=0)
        x_max_data, y_max_data, z_max_data = all_positions.max(axis=0)

        # Ensure minimum range to avoid zero-width axes
        def ensure_min_range(min_val, max_val, min_width=1.0):
            if max_val - min_val < min_width:
                center = (min_val + max_val) / 2
                return center - min_width / 2, center + min_width / 2
            return min_val, max_val

        x_min_data, x_max_data = ensure_min_range(x_min_data, x_max_data)
        y_min_data, y_max_data = ensure_min_range(y_min_data, y_max_data)
        z_min_data, z_max_data = ensure_min_range(z_min_data, z_max_data)

        # Add 10% padding around the data
        padding = 0.1
        x_range = x_max_data - x_min_data
        y_range = y_max_data - y_min_data
        z_range = z_max_data - z_min_data

        x_min = x_min_data - padding * x_range
        x_max = x_max_data + padding * x_range
        y_min = y_min_data - padding * y_range
        y_max = y_max_data + padding * y_range
        z_min = z_min_data - padding * z_range
        z_max = z_max_data + padding * z_range

        # Round to nice numbers for axis labels
        def nice_round(val, direction='down'):
            if abs(val) < 0.001:
                return 0
            magnitude = 10 ** np.floor(np.log10(abs(val)))
            if direction == 'down':
                return np.floor(val / magnitude) * magnitude
            else:
                return np.ceil(val / magnitude) * magnitude

        x_min = nice_round(x_min, 'down')
        x_max = nice_round(x_max, 'up')
        y_min = nice_round(y_min, 'down')
        y_max = nice_round(y_max, 'up')
        z_min = nice_round(z_min, 'down')
        z_max = nice_round(z_max, 'up')

        # Ensure ranges are still valid after rounding
        if x_max <= x_min:
            x_max = x_min + 1
        if y_max <= y_min:
            y_max = y_min + 1
        if z_max <= z_min:
            z_max = z_min + 1

        # Calculate scene scale and zoom based on data extent
        max_extent = max(x_max - x_min, y_max - y_min, z_max - z_min)
        SCENE_SCALE = 10 / max_extent  # Normalize to fit in ~10 units
        VISUAL_SCALE = 1.2  # Scale down axes and grid size
        zoom = 1 / (max_extent * SCENE_SCALE / 10)  # Auto-calculate zoom

        # Camera setup
        self.set_camera_orientation(phi=60 * DEGREES, theta=45 * DEGREES, zoom=zoom)
        self.renderer.camera.light_source.move_to(1000 * IN)
        self.begin_ambient_camera_rotation(rate=0.2)

        # Create axes with actual data coordinates
        axis_step = max_extent / 10  # ~10 tick marks
        axis_step = max(10 ** np.floor(np.log10(axis_step + 0.001)), 0.1)  # Round to power of 10, minimum 0.1

        # Use the same length for all axes (based on max extent)
        uniform_axis_length = max_extent * SCENE_SCALE * VISUAL_SCALE

        # Calculate individual scales for each axis
        x_scale = uniform_axis_length / (x_max_data - x_min_data)
        y_scale = uniform_axis_length / (y_max_data - y_min_data)
        z_scale = uniform_axis_length / (z_max_data - z_min_data)

        axes = ThreeDAxes(
            x_range=[x_min_data, x_max_data, axis_step],
            y_range=[y_min_data, y_max_data, axis_step],
            z_range=[z_min_data, z_max_data, axis_step],
            x_length=uniform_axis_length,
            y_length=uniform_axis_length,
            z_length=uniform_axis_length,
            tips=True,
        )
        self.add(axes)

        # Store bounds for coordinate conversion (use actual data bounds)
        data_bounds = (x_min_data, x_max_data, y_min_data, y_max_data, z_min_data, z_max_data)

        # Create grid aligned with actual data coordinates
        grid = self.create_grid(x_min_data, x_max_data, y_min_data, y_max_data, x_scale, y_scale, axis_step)
        self.add(grid)

        # Set up body visualization
        colors = [RED, GREEN, BLUE, YELLOW, ORANGE, PURPLE, PINK, TEAL]

        stars = []
        paths = []
        path_points = []

        for body in range(NUM_BODIES):
            color = colors[body % len(colors)]
            star = self.create_star(radius=0.15, color=color)

            self.add_fixed_orientation_mobjects(star)

            # Convert data coordinates to scene coordinates
            initial_pos = self.data_to_scene(frames[body][0], data_bounds, x_scale, y_scale, z_scale)
            star.move_to(initial_pos)

            stars.append(star)
            self.add(star)

            # Initialize path points with starting position
            path_points.append([initial_pos.copy()])

            # Create path with initial dummy line (two identical points won't render but avoids errors)
            path = VMobject(stroke_color=color, stroke_width=2, stroke_opacity=0.7, fill_opacity=0)
            path.set_points_as_corners([initial_pos, initial_pos + np.array([0.001, 0, 0])])
            paths.append(path)
            self.add(path)

        # Animation update function
        last_frame = [-1]  # Track last frame to avoid redundant path updates

        def update_all(mob, alpha):
            frame_index = int(alpha * (num_timesteps - 1))

            for body, star in enumerate(stars):
                pos = self.data_to_scene(frames[body][frame_index], data_bounds, x_scale, y_scale, z_scale)
                star.move_to(pos)

            # Only add path points when frame changes
            if frame_index != last_frame[0]:
                last_frame[0] = frame_index
                for body in range(NUM_BODIES):
                    pos = self.data_to_scene(frames[body][frame_index], data_bounds, x_scale, y_scale, z_scale)
                    path_points[body].append(pos.copy())

                    # Update path with at least 2 points
                    if len(path_points[body]) >= 2:
                        paths[body].set_points_as_corners(path_points[body])

        # Run animation - duration based on simulation frames and frame rate
        animation_fps = 24
        run_time = num_timesteps / animation_fps
        
        self.play(
            UpdateFromAlphaFunc(VGroup(*stars), update_all),
            run_time=run_time,
            rate_func=linear
        )

        self.stop_ambient_camera_rotation()
        self.wait(1)

    def data_to_scene(self, pos, data_bounds, x_scale, y_scale, z_scale):
        """Convert simulation coordinates to scene coordinates (centered at origin)."""
        x_min, x_max, y_min, y_max, z_min, z_max = data_bounds
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2
        z_center = (z_min + z_max) / 2
        
        return np.array([
            (pos[0] - x_center) * x_scale,
            (pos[1] - y_center) * y_scale,
            (pos[2] - z_center) * z_scale
        ])

    def create_star(self, radius, color):
        """Create a glowing star effect."""
        # Use Dot instead of Circle to avoid gradient issues in 3D rendering
        core = Dot(radius=radius, color=color, fill_opacity=1)
        halo1 = Dot(radius=radius + 0.02, color=color, fill_opacity=0.25)
        halo2 = Dot(radius=radius + 0.06, color=color, fill_opacity=0.15)
        halo3 = Dot(radius=radius + 0.12, color=color, fill_opacity=0.08)
        return VGroup(halo3, halo2, halo1, core)  # Order: back to front

    def create_grid(self, x_min, x_max, y_min, y_max, x_scale, y_scale, step):
        """Create a grid on the XY plane at z=0."""
        grid = VGroup()
        grid.z_index = -10

        stroke_opacity = 0.2
        stroke_width = 1

        # Calculate center offsets (axes are centered at origin)
        x_center = (x_min + x_max) / 2
        y_center = (y_min + y_max) / 2

        # Find the first tick mark position (multiple of step)
        x_start = np.ceil(x_min / step) * step
        y_start = np.ceil(y_min / step) * step

        # Vertical lines (constant x) - align with axis ticks
        x = x_start
        while x <= x_max + 0.001:
            x_vis = (x - x_center) * x_scale  # Use x_scale
            y_min_vis = (y_min - y_center) * y_scale  # Use y_scale
            y_max_vis = (y_max - y_center) * y_scale  # Use y_scale
            start = np.array([x_vis, y_min_vis, 0])
            end = np.array([x_vis, y_max_vis, 0])
            line = Line(start, end, stroke_color=GREY, stroke_width=stroke_width, stroke_opacity=stroke_opacity)
            grid.add(line)
            x += step

        # Horizontal lines (constant y) - align with axis ticks
        y = y_start
        while y <= y_max + 0.001:
            y_vis = (y - y_center) * y_scale  # Use y_scale
            x_min_vis = (x_min - x_center) * x_scale  # Use x_scale
            x_max_vis = (x_max - x_center) * x_scale  # Use x_scale
            start = np.array([x_min_vis, y_vis, 0])
            end = np.array([x_max_vis, y_vis, 0])
            line = Line(start, end, stroke_color=GREY, stroke_width=stroke_width, stroke_opacity=stroke_opacity)
            grid.add(line)
            y += step

        return grid