import os
from tkinter import IntVar, ttk
import customtkinter
import numpy as np
import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
from pathlib import Path

#leapfrog module
import leapfrog

# Set Appearance of UI
customtkinter.set_appearance_mode("System") 
customtkinter.set_default_color_theme("dark-blue")

#setting cwd directory
CWDDIR = Path.cwd()

class app(customtkinter.CTk):
    def __init__(self):
        customtkinter.CTk.__init__(self)
        
        def get_path(body):
            #Debugging
            if body == -1:
                path = Path(str(CWDDIR)) / "Simulated_Data" #directory
            else:
                path = Path(str(CWDDIR)) / "Simulated_Data" / f"body{body}.csv" #individual csv
            return path

        def show_animation(duration):
            # CONSTANTS
            NUM_BODIES = len(os.listdir(get_path(-1)))
            # length of timeline
            #TIMESTEP = 0.001 to 0.1

            with open(get_path(0), encoding="utf-8") as f:
                row_count = sum(1 for _ in f)

            TIMELINE = np.linspace(0, row_count, row_count)
            # read body positions for every frame from seperate csv files, one for each body
            frames = np.empty((NUM_BODIES, TIMELINE.size, 3), dtype=float)
            for body in range(0, NUM_BODIES):
                path = get_path(body)
                frames[body] = np.genfromtxt(path, delimiter=',')

            # plotting the data
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

            # animation of plot & points
            animated_plots = []
            points = []

            # initialize plots / points on plots
            for body in range(0, len(frames)):
                animated_plots.append(ax.plot([], [], [])[0])
                points.append(ax.plot([], [], [], 'ro', markersize=4)[0])

            def update_data(frame):

                # go through every E.O.M frame by frame (motion has already been calculated)
                for body in range(len(frames)):
                    animated_plots[body].set_data(frames[body, :frame, 0], frames[body, :frame, 1])
                    animated_plots[body].set_3d_properties(frames[body, :frame, 2])

                    points[body].set_data([frames[body, frame, 0]], [frames[body, frame, 1]])
                    points[body].set_3d_properties([frames[body, frame, 2]])

                return animated_plots, points

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


        #logic for custom button
        def custom():
            if var1.get() == 1:
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
            num = int(selection.strip('Stable'))-1
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
            if var1.get()==1:
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

        #submit button logic, checks if variables are valid and enters them in a list
        def update():
            num = int(self.dropdown1.get().strip('Stable'))-1
            precision = 0.01
            if var2.get() == 0:
                precision = round(self.precision.get(), 4)
            elif var2.get() == 1:
                precision = float(self.precisionoverride.get())
            if var1.get() == 1:
                if self.position1.get()[0]!="(" or self.position2.get()[0]!="(" or self.position3.get()[0]!="(" or self.position1.get()[-1]!=")" or self.position2.get()[-1]!=")" or self.position3.get()[-1]!=")" or self.velocity1.get()[0]!="(" or self.velocity2.get()[0]!="(" or self.velocity3.get()[0]!="(" or self.velocity1.get()[-1]!=")" or self.velocity2.get()[-1]!=")" or self.velocity3.get()[-1]!=")":
                    self.textfield.configure(text="Please include AND close brackets for position and velocity")
                    return
                else:
                    customData = list((eval(self.position1.get()),eval(self.mass1.get()),eval(self.velocity1.get()),eval(self.position2.get()),eval(self.mass2.get()),eval(self.velocity2.get()),eval(self.position3.get()),eval(self.mass3.get()),eval(self.velocity3.get())))
                    leapfrog.Simulate(customData, precision, int(self.durationVariable.get()))
                    show_animation(int(self.durationVariable.get()))
            else:
                leapfrog.Simulate(stables[num], precision, int(self.durationVariable.get()))
                show_animation(int(self.durationVariable.get()))

        def override():
            if var2.get() == 1:
                self.precisionoverride.configure(state="normal")
            else:
                self.precisionoverride.configure(state="disabled")

        #stable orbit variable data, currently placeholders
        rt32 = np.sqrt(3)/2
        v = np.sqrt(1/(5*np.sqrt(3)))
        stable1 = list(((5, 0, 0),1,(0,v,0),(-0.5*5,rt32*5,0),1,(-v*rt32,-v/2,0),(-0.5*5,-rt32*5,0),1,(v*rt32,-v/2,0)))
        stable2 = list(((2.57429,0,0),1,(0.216343, 0.332029,0),(-2.57429,0,0),1,(0.216343, 0.332029,0),(0,0,0),1,(-0.432686, -0.664058,0)))
        stable3 = list(((1,0,0),1,(0,0.5, 0),(-0.5,rt32,0),1,(-0.433,-0.25, 0),(-0.5,-rt32,0),1,(0.433,-0.25, 0)))
        stable4 = list(((0,0,20),20,(0,0,10),(10,0,0),3,(20,0,0),(0,0,10),60,(0,0,20)))
        stable5 = list(((1, 3, 0), 3, (0, 0, 0), (-2, -1, 0), 4, (0, 0, 0), (1, -1, 0), 5, (0, 0, 0)))
        stables = list((stable1, stable2, stable3, stable4, stable5))


        #Window config
        self.title("Three Body Problem Simulator")
        self.geometry("1000x800")

        # configure grid layout (4x4)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure((2, 3), weight=0)
        self.grid_rowconfigure((0,3), weight=0)
        self.grid_rowconfigure((1,2), weight=2)

        #drop down for stable systems
        self.dropbox_frame=customtkinter.CTkFrame(self)
        self.dropbox_frame.grid(column=2, row=0, pady=20, padx=20, sticky="nsew")
        self.dropdown1 = customtkinter.CTkOptionMenu(self.dropbox_frame, values=["Stable 1", "Stable 2", "Stable 3", "Stable 4", "Stable 5"], command=stableOrbits)
        self.dropdown1.pack(padx=(20,20), pady=20, side="left")
        #button to allow for custom inputs
        var1 = IntVar()
        self.button = customtkinter.CTkCheckBox(master = self.dropbox_frame, text="Custom", variable=var1, onvalue=1, offvalue=0, corner_radius=8, command=custom)
        self.button.pack(padx=20, pady=20, side="left")

        #tab view for the objects
        self.variables = customtkinter.CTkTabview(self, width=250)
        self.variables.grid(column=2, row=1, rowspan=2, padx=20, pady=(20,0), sticky="nsew")
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
        self.slidertext = customtkinter.CTkLabel(self.submitframe, height=20, width=300, text="")
        self.slidertext.pack(pady=(20,5),side="top")
        self.precision = customtkinter.CTkSlider(self.submitframe, from_=0.1, to=0.0005, width=200, command=lambda value: self.slidertext.configure(text="Precision: " + str(round(value,4)) + " (size of timesteps, lower is more accurate)"), number_of_steps=100)
        self.precision.pack(pady=(5,20), padx=20,side="top")
        self.precision.set(0.01)
        self.slidertext.configure(text="Precision: " + str(round(self.precision.get(),4)) + " (size of timesteps, lower is more accurate)")
        self.durationVariable = customtkinter.CTkEntry(self.submitframe, width=200, placeholder_text="Simulation duration (Seconds)")
        self.durationVariable.pack(pady=20, padx=20,side="top")
        self.submit = customtkinter.CTkButton(self.submitframe,text="Simulate", command=update)
        self.submit.pack(pady=20, padx=20,side="top")
        self.textfield = customtkinter.CTkLabel(self.submitframe, height=20, width=300, text="")
        self.textfield.pack(pady=20,side="top")
        var2= IntVar()
        self.override = customtkinter.CTkCheckBox(self.submitframe, text="Override precision value (this may kill your PC)", variable=var2, onvalue=1, offvalue=0, corner_radius=8, command=override)
        self.override.pack(pady=20, padx=20,side="top")
        self.precisionoverride = customtkinter.CTkEntry(self.submitframe, placeholder_text="Override Precision Value:", width=150)
        self.precisionoverride.configure(state="disabled")
        self.precisionoverride.pack(pady=20, padx=20,side="top")

        #frame for simulation animation
        self.animation_frame = customtkinter.CTkFrame(self, height=500)
        self.animation_frame.grid(column=0, row=0, rowspan=2, columnspan=2, padx=20, pady=20, sticky="nsew")

def main():
    app().mainloop()

if __name__ == "__main__":
    main()
    
# TO DO:
# Add Tests for UI via Pytest
    

