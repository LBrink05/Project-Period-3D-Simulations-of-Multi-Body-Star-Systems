#IMPORT LIBRARIES
import ast
import os
from ast import literal_eval
from tkinter import IntVar, ttk, Toplevel
import customtkinter
import numpy as np
from pathlib import Path
import re 
import datetime
import importlib
import time

#MATPLOTLIB
import matplotlib
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib import colors as mcolors
matplotlib.use("tkagg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation

#CUSTOM MODULES
import data_analysis

#SET GUI THEME
customtkinter.set_appearance_mode("System")
customtkinter.set_default_color_theme("dark-blue")


#SETTING CWD DIRECTORY FOR PATH FUNCTION
CWDDIR = Path.cwd()

#STARTUP TERMINAL MESSAGE
print(" \n--Starting up 3D Multi-Body Stellar Simulator--\n")

#APP
class app(customtkinter.CTk):
    def __init__(self):
        customtkinter.CTk.__init__(self)

        self.popup = None
        self.anim = None
        self.current_fig = None
        self.last_simulation_data = None  # Store last simulation info
        
        # Define all your functions INSIDE __init__ first
        def on_closing():

            # Stop animation
            if self.anim is not None:
                self.anim.event_source.stop()
            
            # Close matplotlib figures
            if self.current_fig is not None:
                plt.close(self.current_fig)
            
            plt.close('all')  # Close any remaining figures
            
            # Close the window

            print("\n--Shutting down 3D Multi-Body Stellar Simulator--\n")
            self.quit()
            self.destroy()
        
        # NOW set the protocol (after the function is defined)
        self.protocol("WM_DELETE_WINDOW", on_closing)

        def get_path(body):
            # Debugging
            if body == -1:
                path = Path(str(CWDDIR)) / "Simulated_Data"  # directory
            else:
                path = Path(str(CWDDIR)) / "Simulated_Data" / f"body{body}.csv"  # individual csv
            return path

        def show_animation(duration):
            # CONSTANTS
            NUM_BODIES = len(os.listdir(get_path(-1)))
            TRAIL=200
            # length of timeline
            #TIMESTEP = 0.001 to 0.1
            with open(get_path(0), encoding="utf-8") as f:
                row_count = sum(1 for _ in f)

            TIMELINE = np.linspace(0, row_count, row_count)
            # read body positions for every frame from seperate csv files, one for each body
            frames = np.empty((NUM_BODIES, TIMELINE.size, 3), dtype=float)
            for body in range(0, NUM_BODIES):
                path = get_path(body)
                frames[body] = np.genfromtxt(path, delimiter=',')[:, :3]

            # plotting the data
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

            #showing time in animation
            plottime = 0
            timetext = ax.text2D(-0.30, 0.955, " t=" + str(plottime), transform=ax.transAxes, fontsize=24, color='black')

            # animation of trails & points (The fading trail was coded using ChatGPT)
            trail_cols = []
            trail_rgbs = []
            points = []

            for body in range(len(frames)):
                # make one fading trail collection per body
                lc = Line3DCollection([], linewidths=2)
                ax.add_collection3d(lc, autolim=False)
                trail_cols.append(lc)

                # store a per-body base color (use default cycle)
                base_color = ax._get_lines.get_next_color()
                trail_rgbs.append(mcolors.to_rgb(base_color))

                # head point (same color as trail)
                points.append(ax.plot([], [], [], 'o', markersize=4, color='red')[0])

            def update_data(frame):
                # go through every E.O.M frame by frame (motion has already been calculated)
                for body in range(len(frames)):
                    start = max(0, frame - TRAIL)
                    xyz = frames[body, start:frame + 1, :]  # (k,3)

                    if xyz.shape[0] < 2:
                        trail_cols[body].set_segments([])
                    else:
                        # segments: (k-1, 2, 3)
                        segs = np.stack([xyz[:-1], xyz[1:]], axis=1)
                        trail_cols[body].set_segments(segs)

                        nseg = segs.shape[0]
                        # alpha from old -> new (0 -> 1)
                        alphas = np.linspace(0.0, 1.0, nseg) ** (3/2)

                        r, g, b = trail_rgbs[body]
                        rgba = np.column_stack([np.full(nseg, r), np.full(nseg, g), np.full(nseg, b), alphas])
                        trail_cols[body].set_color(rgba)

                    # head point
                    points[body].set_data([frames[body, frame, 0]], [frames[body, frame, 1]])
                    points[body].set_3d_properties([frames[body, frame, 2]])

                #getting time passed per frame within simulation
                plottime = float("%.2f"%(frame / 24))
                timetext.set_text("t=" + str(plottime))
                return trail_cols, points, timetext

            self.anim = FuncAnimation(fig, update_data, frames=TIMELINE.size, interval=10, blit=False)

            #destroyes the old animation if ran multiple times
            for widget in self.animation_frame.winfo_children():
                widget.destroy()

            canvas = FigureCanvasTkAgg(fig, self.animation_frame)
            canvas_widget = canvas.get_tk_widget()
            canvas_widget.pack(fill="both", expand=True)

            axis_dim = 25
            ax.set_xlim(-axis_dim, axis_dim)
            ax.set_ylim(-axis_dim, axis_dim)
            ax.set_zlim(-axis_dim, axis_dim)

            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')

            canvas.draw()

        def show_statistics(duration, precision, selection):
            #constants
            NUM_BODIES = len(os.listdir(get_path(-1)))
            
            #opening reference data and simulation data
            path_reference = Path(str(CWDDIR)) / 'Reference_Data' / ((str(selection).split(' - ')[1]) + "_IAS15_dT_" + "0.0005")
            path_simulation = get_path(-1)

             # Check if reference data exists
            if not path_reference.exists():
                print("\n" + "#"*60)
                print("ERROR: Reference data not found!")
                print("#"*60)
                print(f"\nMissing: {path_reference}")
                print("\nPlease run the IAS15 reference simulation first before")
                print("comparing results. Select IAS15 as the integrator and")
                print("run a simulation with precision 0.0005 to generate reference data.")
                return

            # Read phase space data
            reference_data = data_analysis.read_phase_space(NUM_BODIES, path_reference)
            simulated_data = data_analysis.read_phase_space(NUM_BODIES, path_simulation)
            
            # Extract masses from configuration
            if selection.startswith('S'):
                num = int(selection.split('-')[0].strip().split(' ')[1]) - 1
                masses = [stables[num][1], stables[num][4], stables[num][7]]
            else:
                num = int(selection.split('-')[0].strip().split(' ')[1]) - 1
                masses = [customs[num][1], customs[num][4], customs[num][7]]
            
            # Create error function
            error_func = data_analysis.error_function(reference_data, simulated_data, masses)

            # count frames - use minimum length between reference and simulation
            min_frames = min(len(simulated_data[0]), len(reference_data[0]))

            # timeline: 1 time unit = 24 frames
            frames = np.arange(min_frames)
            TIMELINE = frames / 24

            # Calculate errors at each time point
            trajectory_errors = []
            hamiltonian_errors = []
            
            for frame in frames:
                traj_err, ham_err = error_func(frame)
                trajectory_errors.append(traj_err)
                hamiltonian_errors.append(ham_err)
            
            # Calculate max errors
            traj_max = data_analysis.calculate_max_error(trajectory_errors, dt=1.0/24.0)
            ham_max = data_analysis.calculate_max_error(hamiltonian_errors, dt=1.0/24.0)

            # Create plot with two subplots
            fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
            
            # Trajectory Error Plot
            ax1.plot(TIMELINE, trajectory_errors, linewidth=2, color='red', label='Trajectory Error')
            ax1.axhline(y=traj_max, color='darkred', linestyle='--', linewidth=1.5, 
                       label=f'%_MAX = {traj_max:.2e}%')
            ax1.set_xlabel('Time', fontsize=11)
            ax1.set_ylabel('Trajectory Error E_% (%)', fontsize=11)
            ax1.set_title('Phase Space Trajectory Error vs Time', fontsize=12, fontweight='bold')
            ax1.grid(True, alpha=0.3, which='both')
            ax1.legend(loc='upper right')
            
            # Hamiltonian Error Plot
            ax2.plot(TIMELINE, hamiltonian_errors, linewidth=2, color='blue', label='Hamiltonian Error')
            ax2.axhline(y=ham_max, color='darkblue', linestyle='--', linewidth=1.5,
                       label=f'H_MAX = {ham_max:.2e}%')
            ax2.set_xlabel('Time', fontsize=11)
            ax2.set_ylabel('Hamiltonian Error E_H% (%)', fontsize=11)
            ax2.set_title('Energy Conservation Error vs Time', fontsize=12, fontweight='bold')
            ax2.grid(True, alpha=0.3, which='both')
            ax2.legend(loc='upper right')
            
            plt.tight_layout()
            dt_str = f"{precision:g}".replace(".", "p")  # 0.065 -> "0p065"
            out = Path(str(CWDDIR)) / "Statistics" / f"{self.dropdown1.get()}_Simulated_data_dT_{dt_str}.png"
            out.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(out, dpi=300, bbox_inches="tight")

            # Destroy old statistics plot if it exists
            for widget in self.statistics_frame.winfo_children():
                widget.destroy()
            
            # Embed in GUI
            canvas2 = FigureCanvasTkAgg(fig2, self.statistics_frame)
            canvas_widget2 = canvas2.get_tk_widget()
            canvas_widget2.pack(fill="both", expand=True)
            canvas2.draw()

            # Print summary statistics to console
            print("\n" + "#"*60)
            print(f"ERROR ANALYSIS: ")
            print("#"*60 + "\n")
            print(f"Total frames analyzed: {min_frames}")
            print(f"Max Trajectory Error: {traj_max:.4e}%")
            print(f"Max Hamiltonian Error: {ham_max:.4e}%\n")

        def show_lyapunov_analysis():
            
            if self.last_simulation_data is None:
                self.lyapunov_status.configure(text="Please run a simulation first")
                return
            
            # Update status
            self.lyapunov_status.configure(text="Calculating Lyapunov exponents...")
            self.update()
            
            # Extract simulation parameters
            NUM_BODIES = self.last_simulation_data['num_bodies']
            masses = self.last_simulation_data['masses']
            precision = self.last_simulation_data['precision']
            
            # Read simulated data
            path_simulation = get_path(-1)
            simulated_data = data_analysis.read_phase_space(NUM_BODIES, path_simulation)
            
            # Calculate Lyapunov exponents
            start_time = time.time()
            lyapunov_spectrum, lyapunov_time, exponents_over_time, time_points, sorted_indices = \
                data_analysis.calculate_lyapunov_exponents(simulated_data, masses, dt=precision, renorm_interval=10)
            calc_time = time.time() - start_time
            
            self.lyapunov_status.pack_forget()

            # Create visualization
            fig3 = plt.figure(figsize=(12, 10))
            
            # Plot 1: Full Lyapunov Spectrum (bar plot)
            ax1 = plt.subplot(3, 1, 1)
            indices = np.arange(len(lyapunov_spectrum))
            colors_spectrum = ['red' if x > 0 else 'blue' if x < 0 else 'gray' for x in lyapunov_spectrum]
            ax1.bar(indices, lyapunov_spectrum, color=colors_spectrum, alpha=0.7, edgecolor='black')
            ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
            
            # Create dimension labels for x-axis based on original (unsorted) indices
            # For 3 bodies: x1,y1,z1,vx1,vy1,vz1, x2,y2,z2,vx2,vy2,vz2, x3,y3,z3,vx3,vy3,vz3
            all_dimension_labels = []
            coord_names = ['x', 'y', 'z', 'vx', 'vy', 'vz']
            for body in range(1, NUM_BODIES + 1):
                for coord in coord_names:
                    all_dimension_labels.append(f'{coord}{body}')
            
            # Map sorted indices to their original dimension labels
            dimension_labels_sorted = [all_dimension_labels[idx] for idx in sorted_indices]
            
            # Set x-axis with dimension labels
            ax1.set_xticks(indices)
            ax1.set_xticklabels(dimension_labels_sorted, rotation=45, ha='right', fontsize=8)
            ax1.set_xlabel('Phase Space Dimension (sorted by λ value)', fontsize=11)
            ax1.set_ylabel('Lyapunov Exponent λ_i', fontsize=11)
            ax1.set_title('Full Lyapunov Spectrum (sorted descending)', fontsize=12, fontweight='bold')
            ax1.grid(True, alpha=0.3, axis='y')
            
            # Add text annotation for max exponent
            if len(lyapunov_spectrum) > 0:
                ax1.text(0.02, 0.98, f'λ_max = {lyapunov_spectrum[0]:.4e} ({dimension_labels_sorted[0]})\nt_L = {lyapunov_time:.4e}',
                        transform=ax1.transAxes, fontsize=10, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            # Plot 2: Evolution of top 6 exponents over renormalization steps
            ax2 = plt.subplot(3, 1, 2)
            if len(exponents_over_time) > 0 and len(time_points) > 0:
                exponents_array = np.array(exponents_over_time)
                # Get indices of top 3 (most positive) and bottom 3 (most negative)
                top_3_indices = list(range(3))
                bottom_3_indices = list(range(len(lyapunov_spectrum) - 3, len(lyapunov_spectrum)))
                indices_to_plot = top_3_indices + bottom_3_indices
                # Distinct colors for each exponent
                exponent_colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628']
                for idx, i in enumerate(indices_to_plot):
                    ax2.scatter(time_points, exponents_array[:, i], 
                                label=f'{dimension_labels_sorted[i]}: {lyapunov_spectrum[i]:.3e}',
                                s=50, color=exponent_colors[idx], alpha=0.8)
                ax2.axhline(y=0, color='black', linestyle='--', linewidth=0.8)
                ax2.set_xlabel('Renormalization Steps', fontsize=11)
                ax2.set_ylabel('Lyapunov Exponent', fontsize=11)
                ax2.set_title('Evolution of Top 6 Lyapunov Exponents (by magnitude)', fontsize=12, fontweight='bold')
                ax2.legend(loc='best', fontsize=9)
                ax2.grid(True, alpha=0.3)

            # Plot 3: Phase space volume preservation (sum of exponents)
            ax3 = plt.subplot(3, 1, 3)
            if len(exponents_over_time) > 0:
                exponents_array = np.array(exponents_over_time)
                sum_exponents = np.sum(exponents_array, axis=1)
                ax3.scatter(time_points, sum_exponents, s=50, color='green', label='Σλ_i')
                ax3.axhline(y=0, color='black', linestyle='--', linewidth=0.8)
                ax3.set_xlabel('Renormalization Steps', fontsize=11)
                ax3.set_ylabel('Sum of Exponents', fontsize=11)
                ax3.set_title('Phase Space Volume Preservation (Liouville\'s Theorem)', fontsize=12, fontweight='bold')
                ax3.grid(True, alpha=0.3)
                # Annotation about Liouville's theorem
                final_sum = sum_exponents[-1] if len(sum_exponents) > 0 else 0
                ax3.text(0.02, 0.98, f'Final Σλ_i = {final_sum:.4e}\n(Should be ≈ 0 for Hamiltonian systems)',
            transform=ax3.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
            
            plt.tight_layout()
            
            # Save figure
            dt_str = f"{precision:g}".replace(".", "p")
            out = Path(str(CWDDIR)) / "Statistics" / f"{self.dropdown1.get()}_Lyapunov_dT_{dt_str}.png"
            out.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(out, dpi=300, bbox_inches="tight")
            
            # Clear old plot and display new one
            for widget in self.lyapunov_frame.winfo_children():
                widget.destroy()
            
            canvas3 = FigureCanvasTkAgg(fig3, self.lyapunov_frame)
            canvas_widget3 = canvas3.get_tk_widget()
            canvas_widget3.pack(fill="both", expand=True)
            canvas3.draw()
            
            # Print detailed results to console
            print("\n" + "#"*60)
            print("LYAPUNOV EXPONENT ANALYSIS")
            print("#"*60 + "\n")
            print(f"Configuration: {self.dropdown1.get()}")
            print(f"Timestep (dt): {precision}")
            print(f"Calculation time: {calc_time:.2f} seconds")
            print(f"\nMaximum Lyapunov Exponent: {lyapunov_spectrum[0]:.6e}")
            print(f"Lyapunov Time (t_L = 1/λ_max): {lyapunov_time:.6e}")
            
            if lyapunov_spectrum[0] > 1e-6:
                print(f"\nSystem is CHAOTIC (λ_max > 0)")
            elif lyapunov_spectrum[0] > -1e-6:
                print(f"\nSystem is MARGINALLY STABLE (λ_max ≈ 0)")
            else:
                print(f"\nSystem is STABLE (λ_max < 0)")
            
            print(f"\nNumber of positive exponents: {np.sum(lyapunov_spectrum > 0)}")
            print(f"Number of negative exponents: {np.sum(lyapunov_spectrum < 0)}")
            print(f"Number of near-zero exponents: {np.sum(np.abs(lyapunov_spectrum) < 1e-6)}")
            
            if len(lyapunov_spectrum) > 0:
                sum_lambda = np.sum(lyapunov_spectrum)
                print(f"\nSum of all exponents (Liouville): {sum_lambda:.6e}")
                print(f"  (Should be ≈ 0 for Hamiltonian systems)")

        #logic for custom button
        def custom():
            if self.var1.get() == 1:
                self.position1.configure(state="normal")
                self.mass1.configure(state="normal")
                self.velocity1.configure(state="normal")
                self.position2.configure(state="normal")
                self.mass2.configure(state="normal")
                self.velocity2.configure(state="normal")
                self.position3.configure(state="normal")
                self.mass3.configure(state="normal")
                self.velocity3.configure(state="normal")
            else:
                #resets variables to selected stable orbit
                stableOrbits(self.dropdown1.get())

                self.position1.configure(state="disabled")
                self.mass1.configure(state="disabled")
                self.velocity1.configure(state="disabled")
                self.position2.configure(state="disabled")
                self.mass2.configure(state="disabled")
                self.velocity2.configure(state="disabled")
                self.position3.configure(state="disabled")
                self.mass3.configure(state="disabled")
                self.velocity3.configure(state="disabled")

        #logic for dropdown selection and value filling
        def stableOrbits(selection):
            if selection.startswith('C'):
                num = int(selection.split('-')[0].strip().split(' ')[1]) - 1
                self.position1.configure(state="normal")
                self.position1.delete(0, "end")
                self.position1.insert(0, "(" + str(round(customs[num][0][0], 3)) + "," + str(round(customs[num][0][1], 3)) + "," + str(round(customs[num][0][2], 3)) + ")")

                self.mass1.configure(state="normal")
                self.mass1.delete(0, "end")
                self.mass1.insert(0, customs[num][1])

                self.velocity1.configure(state="normal")
                self.velocity1.delete(0, "end")
                self.velocity1.insert(0, "(" + str(round(customs[num][2][0], 3)) + "," + str(round(customs[num][2][1], 3)) + "," + str(round(customs[num][2][2], 3)) + ")")

                self.position2.configure(state="normal")
                self.position2.delete(0, "end")
                self.position2.insert(0, "(" + str(round(customs[num][3][0], 3)) + "," + str(round(customs[num][3][1], 3)) + "," + str(round(customs[num][3][2], 3)) + ")")

                self.mass2.configure(state="normal")
                self.mass2.delete(0, "end")
                self.mass2.insert(0, customs[num][4])

                self.velocity2.configure(state="normal")
                self.velocity2.delete(0, "end")
                self.velocity2.insert(0, "(" + str(round(customs[num][5][0], 3)) + "," + str(round(customs[num][5][1], 3)) + "," + str(round(customs[num][5][2], 3)) + ")")

                self.position3.configure(state="normal")
                self.position3.delete(0, "end")
                self.position3.insert(0, "(" + str(round(customs[num][6][0], 3)) + "," + str(round(customs[num][6][1], 3)) + "," + str(round(customs[num][6][2], 3)) + ")")

                self.mass3.configure(state="normal")
                self.mass3.delete(0, "end")
                self.mass3.insert(0, customs[num][7])

                self.velocity3.configure(state="normal")
                self.velocity3.delete(0, "end")
                self.velocity3.insert(0, "(" + str(round(customs[num][8][0], 3)) + "," + str(round(customs[num][8][1], 3)) + "," + str(round(customs[num][8][2], 3)) + ")")
            else:
                num = int(selection.split('-')[0].strip().split(' ')[1]) - 1
                self.position1.configure(state="normal")
                self.position1.delete(0,"end")
                self.position1.insert(0, "("+str(round(stables[num][0][0],3))+","+str(round(stables[num][0][1],3))+","+str(round(stables[num][0][2],3))+")")

                self.mass1.configure(state="normal")
                self.mass1.delete(0,"end")
                self.mass1.insert(0, stables[num][1])

                self.velocity1.configure(state="normal")
                self.velocity1.delete(0,"end")
                self.velocity1.insert(0, "("+str(round(stables[num][2][0],3))+","+str(round(stables[num][2][1],3))+","+str(round(stables[num][2][2],3))+")")

                self.position2.configure(state="normal")
                self.position2.delete(0,"end")
                self.position2.insert(0, "("+str(round(stables[num][3][0],3))+","+str(round(stables[num][3][1],3))+","+str(round(stables[num][3][2],3))+")")

                self.mass2.configure(state="normal")
                self.mass2.delete(0,"end")
                self.mass2.insert(0, stables[num][4])

                self.velocity2.configure(state="normal")
                self.velocity2.delete(0,"end")
                self.velocity2.insert(0, "("+str(round(stables[num][5][0],3))+","+str(round(stables[num][5][1],3))+","+str(round(stables[num][5][2],3))+")")

                self.position3.configure(state="normal")
                self.position3.delete(0,"end")
                self.position3.insert(0, "("+str(round(stables[num][6][0],3))+","+str(round(stables[num][6][1],3))+","+str(round(stables[num][6][2],3))+")")

                self.mass3.configure(state="normal")
                self.mass3.delete(0,"end")
                self.mass3.insert(0, stables[num][7])

                self.velocity3.configure(state="normal")
                self.velocity3.delete(0,"end")
                self.velocity3.insert(0, "("+str(round(stables[num][8][0],3))+","+str(round(stables[num][8][1],3))+","+str(round(stables[num][8][2],3))+")")

                #keeps variables editable/disabled based on checkbox status
            if self.var1.get()==1:
                return
            else:
                self.position1.configure(state="disabled")
                self.mass1.configure(state="disabled")
                self.velocity1.configure(state="disabled")
                self.position2.configure(state="disabled")
                self.mass2.configure(state="disabled")
                self.velocity2.configure(state="disabled")
                self.position3.configure(state="disabled")
                self.mass3.configure(state="disabled")
                self.velocity3.configure(state="disabled")

        def integrator(selection):
            self.dropdown2.set(selection)
            
        #submits button logic, checks if variables are valid and enters them in a list
        def update():
            
            #get time when button was pressed
            submission_time = datetime.datetime.now()

            selection = self.dropdown1.get()
            integrator_selection = importlib.import_module(self.dropdown2.get())

            print("Using Numerical Method: ", str(integrator_selection) + "\n")
            
            num = int(selection.split('-')[0].strip().split(' ')[1]) - 1
            precision = 0.01
            if self.var2.get() == 0:
                precision = round(self.precision.get(), 4)
            elif self.var2.get() == 1:
                precision = float(self.precisionoverride.get())
            
            # Extract masses for storage
            if selection.startswith('S'):
                masses = [stables[num][1], stables[num][4], stables[num][7]]
            else:
                masses = [customs[num][1], customs[num][4], customs[num][7]]
            
            if self.var1.get() == 1:
                if self.position1.get()[0]!="(" or self.position2.get()[0]!="(" or self.position3.get()[0]!="(" or self.position1.get()[-1]!=")" or self.position2.get()[-1]!=")" or self.position3.get()[-1]!=")" or self.velocity1.get()[0]!="(" or self.velocity2.get()[0]!="(" or self.velocity3.get()[0]!="(" or self.velocity1.get()[-1]!=")" or self.velocity2.get()[-1]!=")" or self.velocity3.get()[-1]!=")":
                    self.textfield.configure(text="Please include AND close brackets for position and velocity")
                    return
                else:
                    customData = list((eval(self.position1.get()),eval(self.mass1.get()),eval(self.velocity1.get()),eval(self.position2.get()),eval(self.mass2.get()),eval(self.velocity2.get()),eval(self.position3.get()),eval(self.mass3.get()),eval(self.velocity3.get())))
                    masses = [customData[1], customData[4], customData[7]]
                    integrator_selection.Simulate(customData, precision, int(self.durationVariable.get()))
                    show_animation(int(self.durationVariable.get()))
            else:
                if self.dropdown1.get().startswith('S'):
                    integrator_selection.Simulate(stables[num], precision, int(self.durationVariable.get()))
                    show_animation(int(self.durationVariable.get()))
                else:
                    integrator_selection.Simulate(customs[num], precision, int(self.durationVariable.get()))
                    show_animation(int(self.durationVariable.get()))
            
            rendering_time = datetime.datetime.now()
            time_elapsed = rendering_time -  submission_time
            print("Processing time of simulation: " + str(time_elapsed) + "\n")
            
            # Store simulation data for Lyapunov analysis
            NUM_BODIES = 3
            self.last_simulation_data = {
                'num_bodies': NUM_BODIES,
                'masses': masses,
                'precision': precision,
                'duration': int(self.durationVariable.get()),
                'selection': selection
            }

            show_statistics(int(self.durationVariable.get()), precision, selection)

        def override():
            if self.var2.get() == 1:
                self.precisionoverride.configure(state="normal")
            else:
                self.precisionoverride.configure(state="disabled")

        def save():
            if self.var1.get()!=1:
                return

            def submit():
                new = list((ast.literal_eval(self.position1.get()), eval(self.mass1.get()), eval(self.velocity1.get()),
                            eval(self.position2.get()), eval(self.mass2.get()), eval(self.velocity2.get()),
                            eval(self.position3.get()), eval(self.mass3.get()), eval(self.velocity3.get()), self.popup.name.get().split('- ')[-1]))
                if self.popup.new_save.get():
                    presets = open("Configurations/Customs.csv", "a")
                    customs.append(new)
                    presets.write(str(new) + '\n')

                    full.append(f'Custom {len(full) - len(stables) + 1} - {new[-1]}')
                    self.dropdown1.set(full[-1])
                    presets.close()
                else:
                    if self.dropdown1.get().startswith("S"):
                        print("Can not edit stable orbits")
                        return
                    index = int(self.dropdown1.get().split('-')[0].strip().split(' ')[1]) - 1
                    customs[index] = new
                    presets = open("Configurations/Customs.csv", "w")
                    presets.write('\n')
                    for i in customs:
                        presets.write(str(i) + '\n')

                    full[len(stables) + index] = f'Custom {index + 1} - {new[-1]}'
                    self.dropdown1.set(full[len(stables) + index])

                self.dropdown1.configure(values=full)

                self.popup.destroy()
                self.popup.update()
                return

            self.popup = customtkinter.CTkToplevel(self)
            self.popup.title("Name Config")
            self.popup.geometry("400x200")
            self.popup.resizable(False, False)
            self.popup.attributes("-topmost", True)

            self.popup.grid_columnconfigure(1, weight=1)
            self.popup.grid_rowconfigure(1, weight=1)

            self.popup.frame = customtkinter.CTkFrame(self.popup)
            self.popup.frame.grid(row=1, column=1, rowspan=2, columnspan=2, pady=0, padx=0, sticky="nsew")
            self.popup.new_save = customtkinter.CTkCheckBox(self.popup.frame, text="New Save", onvalue=1, offvalue=0, command=check)
            self.popup.new_save.pack(side="left", padx=(40,5), pady=5)
            self.popup.name = customtkinter.CTkEntry(self.popup.frame, placeholder_text="Save Name")
            self.popup.name.insert(0,self.dropdown1.get())
            self.popup.name.configure(state="disabled")
            self.popup.name.pack(side="right", padx=(5,40), pady=5)

            self.popup.submit = customtkinter.CTkButton(self.popup.frame, text="Submit", command=submit)
            self.popup.submit.pack(side="bottom", padx=5, pady=(5,20))

        def check():
            if self.popup.new_save.get():
                self.popup.name.configure(state="normal")

                self.popup.name.delete(0,"end")
            else:
                self.popup.name.delete(0,"end")
                self.popup.name.insert(0,self.dropdown1.get())
                self.popup.name.configure(state="disabled")
        
        # lists of configurations
        stables = list()
        customs = list()
        full = list()

        # lists of numerical methods
        integrator_list = list(["leapfrog", "forwardeuler", "IAS15", "Adaptive_leapfrog"])

        with open("Configurations/Stables.csv", "r") as presets:
            for x in presets:
                line = x.strip()
                if line:
                    stables.append(ast.literal_eval(line))

        with open("Configurations/Customs.csv", "r") as presets:
            for x in presets:
                line = x.strip()
                if line:
                    customs.append(ast.literal_eval(line))


        for i, stable in enumerate(stables):
            full.append(f'Stable {i+1} - {stable[-1]}') 
        for x, customized in enumerate(customs):
            full.append(f'Custom {x+1} - {customized[-1]}')

        #Window config
        self.title("Three Body Problem Simulator")
        self.geometry("1400x1000")

        # configure grid layout (4x4)
        self.grid_columnconfigure((1), weight=1)
        self.grid_columnconfigure((2, 3), weight=0)
        self.grid_rowconfigure(0, weight=0)
        self.grid_rowconfigure((1, 2, 3), weight=1)

        #tab views for simulations
        self.simtabs = customtkinter.CTkTabview(self, width=300, height=800)
        self.simtabs.grid(column=1, row=0, rowspan=4, columnspan=1, padx=20, pady=20, sticky="nsew")
        self.simtabs.add("Simulation")
        self.simtabs.add("Lyapunov")
        #self.simtabs.add("tab 3")
        self.simtabs.tab("Simulation").grid_columnconfigure((0,1), weight=1)
        self.simtabs.tab("Simulation").grid_rowconfigure((0,2), weight=1)
        self.simtabs.tab("Lyapunov").grid_columnconfigure((0,1), weight=1)
        self.simtabs.tab("Lyapunov").grid_rowconfigure(0, weight=0)  # Control frame - fixed height
        self.simtabs.tab("Lyapunov").grid_rowconfigure(1, weight=1)  # Plot frame - expandable
        #self.simtabs.tab("tab 3").grid_columnconfigure((0,1), weight=1)
        #self.simtabs.tab("tab 3").grid_rowconfigure((0,1,2), weight=1)
        
        #drop down for stable systems
        self.dropbox_frame=customtkinter.CTkFrame(self, height=300)
        self.dropbox_frame.grid(column=2,row=0, pady=(20,5), padx=20, sticky="nsew")
        self.dropdown1 = customtkinter.CTkOptionMenu(self.dropbox_frame, values=full, command=stableOrbits)
        self.dropdown1.pack(padx=(20,20), pady=20, side="left")

        #drop down for numerical integrators
        self.dropdown2 = customtkinter.CTkOptionMenu(self.dropbox_frame, values=integrator_list, command=integrator)
        self.dropdown2.pack(padx=(20,20), pady=20, side="left")

        #button to allow for custom inputs
        self.var1 = IntVar()
        self.button = customtkinter.CTkCheckBox(master = self.dropbox_frame, text="Custom", variable=self.var1, onvalue=1, offvalue=0, corner_radius=8, command=custom)
        self.button.pack(padx=20, pady=20, side="left")

        #tab view for the objects
        self.variables = customtkinter.CTkTabview(self, width=250, height=400)
        self.variables.grid(column=2, row=1, rowspan=2, padx=20, pady=(10,0), sticky="nsew")
        self.variables.add("Object 1")
        self.variables.add("Object 2")
        self.variables.add("Object 3")
        self.variables.tab("Object 1").grid_columnconfigure((0,1), weight=1)
        self.variables.tab("Object 1").grid_rowconfigure((0,1,2), weight=1)
        self.variables.tab("Object 2").grid_columnconfigure((0,1), weight=1)
        self.variables.tab("Object 2").grid_rowconfigure((0,1,2), weight=1)
        self.variables.tab("Object 3").grid_columnconfigure((0,1), weight=1)
        self.variables.tab("Object 3").grid_rowconfigure((0,1,2), weight=1)

        #Object 1 variables
        self.positionLabel1 = customtkinter.CTkLabel(self.variables.tab("Object 1"), text="Position (X,Y,Z):")
        self.positionLabel1.grid(column=0, row=0, padx=20, pady=20, sticky="nsew")
        self.position1 = customtkinter.CTkEntry(self.variables.tab("Object 1"), state="disabled")
        self.position1.grid(column=1, row=0, padx=20, pady=20)
        self.massLabel1 = customtkinter.CTkLabel(self.variables.tab("Object 1"), text="Mass (Solar Masses):")
        self.massLabel1.grid(column=0, row=1, padx=20, pady=20, sticky="nsew")
        self.mass1 = customtkinter.CTkEntry(self.variables.tab("Object 1"), state="disabled")
        self.mass1.grid(column=1, row=1, padx=20, pady=20)
        self.velocityLabel1 = customtkinter.CTkLabel(self.variables.tab("Object 1"), text="Velocity (km/s) (X,Y,Z):")
        self.velocityLabel1.grid(column=0, row=2, padx=20, pady=20, sticky="nsew")
        self.velocity1 = customtkinter.CTkEntry(self.variables.tab("Object 1"), state="disabled")
        self.velocity1.grid(column=1, row=2, padx=20, pady=20)

        #Object 2 variables
        self.positionLabel2 = customtkinter.CTkLabel(self.variables.tab("Object 2"), text="Position (X,Y,Z):")
        self.positionLabel2.grid(column=0, row=0, padx=20, pady=20, sticky="nsew")
        self.position2 = customtkinter.CTkEntry(self.variables.tab("Object 2"), state="disabled")
        self.position2.grid(column=1, row=0, padx=20, pady=20)
        self.massLabel2 = customtkinter.CTkLabel(self.variables.tab("Object 2"), text="Mass (Solar Masses):")
        self.massLabel2.grid(column=0, row=1, padx=20, pady=20, sticky="nsew")
        self.mass2 = customtkinter.CTkEntry(self.variables.tab("Object 2"), state="disabled")
        self.mass2.grid(column=1, row=1, padx=20, pady=20)
        self.velocityLabel2 = customtkinter.CTkLabel(self.variables.tab("Object 2"), text="Velocity (km/s) (X,Y,Z):")
        self.velocityLabel2.grid(column=0, row=2, padx=20, pady=20, sticky="nsew")
        self.velocity2 = customtkinter.CTkEntry(self.variables.tab("Object 2"), state="disabled")
        self.velocity2.grid(column=1, row=2, padx=20, pady=20)

        # Object 3 variables
        self.positionLabel3 = customtkinter.CTkLabel(self.variables.tab("Object 3"), text="Position (X,Y,Z):")
        self.positionLabel3.grid(column=0, row=0, padx=20, pady=20, sticky="nsew")
        self.position3 = customtkinter.CTkEntry(self.variables.tab("Object 3"), state="disabled")
        self.position3.grid(column=1, row=0, padx=20, pady=20)
        self.massLabel3 = customtkinter.CTkLabel(self.variables.tab("Object 3"), text="Mass (Solar Masses):")
        self.massLabel3.grid(column=0, row=1, padx=20, pady=20, sticky="nsew")
        self.mass3 = customtkinter.CTkEntry(self.variables.tab("Object 3"), state="disabled")
        self.mass3.grid(column=1, row=1, padx=20, pady=20)
        self.velocityLabel3 = customtkinter.CTkLabel(self.variables.tab("Object 3"), text="Velocity (km/s) (X,Y,Z):")
        self.velocityLabel3.grid(column=0, row=2, padx=20, pady=20, sticky="nsew")
        self.velocity3 = customtkinter.CTkEntry(self.variables.tab("Object 3"), state="disabled")
        self.velocity3.grid(column=1, row=2, padx=20, pady=20)

        #Submit button and hidden textfield for error popups
        self.submitframe = customtkinter.CTkFrame(self)
        self.submitframe.grid(column=2, row=3, pady=(0,20), padx=20, sticky="nsew")
        self.save = customtkinter.CTkButton(self.submitframe, height=20, width=150, text="Save Configuration", command=save)
        self.save.pack(pady=(5,5), side="top")
        self.slidertext = customtkinter.CTkLabel(self.submitframe, height=20, width=300, text="")
        self.slidertext.pack(pady=(20,5),side="top")
        self.precision = customtkinter.CTkSlider(self.submitframe, from_=0.1, to=0.0005, width=200, command=lambda value: self.slidertext.configure(text="Precision: " + str(round(value,4)) + " (size of timesteps, lower is more accurate)"), number_of_steps=100)
        self.precision.pack(pady=(5,20), padx=20,side="top")
        self.precision.set(0.01)
        self.slidertext.configure(text="Precision: " + str(round(self.precision.get(),4)) + " (size of timesteps, lower is more accurate)")
        self.durationVariable = customtkinter.CTkEntry(self.submitframe, width=200, placeholder_text="Simulation duration (s)")
        self.durationVariable.pack(pady=10, padx=20,side="top")
        self.submit = customtkinter.CTkButton(self.submitframe,text="Simulate", command=update)
        self.submit.pack(pady=10, padx=20,side="top")
        self.textfield = customtkinter.CTkLabel(self.submitframe, height=20, width=300, text="")
        self.textfield.pack(pady=10,side="top")
        self.var2= IntVar()
        self.override = customtkinter.CTkCheckBox(self.submitframe, text="Override precision value (this may kill your PC)", variable=self.var2, onvalue=1, offvalue=0, corner_radius=8, command=override)
        self.override.pack(pady=20, padx=20,side="top")
        self.precisionoverride = customtkinter.CTkEntry(self.submitframe, placeholder_text="Override Precision Value:", width=150)
        self.precisionoverride.configure(state="disabled")
        self.precisionoverride.pack(pady=20, padx=20,side="top")

        #frame for simulation animation
        self.animation_frame = customtkinter.CTkFrame(self.simtabs.tab("Simulation"), height=500)
        self.animation_frame.grid(column=0, row=0, rowspan=2, columnspan=2, padx=20, pady=20, sticky="nsew")

        #frame for statistics plot
        self.statistics_frame = customtkinter.CTkFrame(self.simtabs.tab("Simulation"), height=300)
        self.statistics_frame.grid(column=0, row=2, rowspan=2, columnspan=2, padx=20, pady=(0,20), sticky="nsew")

        # Lyapunov tab setup - control frame with fixed height
        self.lyapunov_control_frame = customtkinter.CTkFrame(self.simtabs.tab("Lyapunov"), height=100)
        self.lyapunov_control_frame.grid(column=0, row=0, columnspan=2, padx=20, pady=(20,10), sticky="ew")
        self.lyapunov_control_frame.grid_propagate(False)  # Prevent frame from shrinking
        
        self.lyapunov_button = customtkinter.CTkButton(self.lyapunov_control_frame, text="Calculate Lyapunov Exponents", 
                                                       command=show_lyapunov_analysis, height=40, font=("Arial", 14, "bold"))
        self.lyapunov_button.pack(pady=15, padx=20)
        
        self.lyapunov_status = customtkinter.CTkLabel(self.lyapunov_control_frame, text="Run a simulation first, then click to analyze", 
                                                      font=("Arial", 11))
        self.lyapunov_status.pack(pady=(0,10))
        
        # Lyapunov plot frame - expandable
        self.lyapunov_frame = customtkinter.CTkFrame(self.simtabs.tab("Lyapunov"), height=600)
        self.lyapunov_frame.grid(column=0, row=1, columnspan=2, padx=20, pady=(0,20), sticky="nsew")

        integrator("leapfrog")
        stableOrbits("Stable 1 - Equilateral Triangle")

app().mainloop()