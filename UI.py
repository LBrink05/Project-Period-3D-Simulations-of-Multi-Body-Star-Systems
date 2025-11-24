from tkinter import IntVar
import customtkinter
import tkinter
import tkinter.messagebox as messagebox

from customtkinter import CTk, CTkEntry

customtkinter.set_appearance_mode("System")
customtkinter.set_default_color_theme("dark-blue")

class app(customtkinter.CTk):
    def __init__(self):
        customtkinter.CTk.__init__(self)

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
                self.position1.configure(state="disabled")
                self.mass1.configure(state="disabled")
                self.velocity1.configure(state="disabled")
                self.position2.configure(state="disabled")
                self.mass2.configure(state="disabled")
                self.velocity2.configure(state="disabled")
                self.position3.configure(state="disabled")
                self.mass3.configure(state="disabled")
                self.velocity3.configure(state="disabled")

        def stableOrbits(selection):
            num = int(selection.strip('Stable'))
            num = num-1
            self.position1.configure(state="normal")
            self.position1.delete(0,"end")
            self.position1.insert(0, "("+str(stables[num][0][0])+","+str(stables[num][0][1])+","+str(stables[num][0][2])+")")

            self.mass1.configure(state="normal")
            self.mass1.delete(0,"end")
            self.mass1.insert(0, stables[num][1])

            self.velocity1.configure(state="normal")
            self.velocity1.delete(0,"end")
            self.velocity1.insert(0, "("+str(stables[num][2][0])+","+str(stables[num][2][1])+","+str(stables[num][2][2])+")")

            self.position2.configure(state="normal")
            self.position2.delete(0,"end")
            self.position2.insert(0, "("+str(stables[num][3][0])+","+str(stables[num][3][1])+","+str(stables[num][3][2])+")")

            self.mass2.configure(state="normal")
            self.mass2.delete(0,"end")
            self.mass2.insert(0, stables[num][4])

            self.velocity2.configure(state="normal")
            self.velocity2.delete(0,"end")
            self.velocity2.insert(0, "("+str(stables[num][5][0])+","+str(stables[num][5][1])+","+str(stables[num][5][2])+")")

            self.position3.configure(state="normal")
            self.position3.delete(0,"end")
            self.position3.insert(0, "("+str(stables[num][6][0])+","+str(stables[num][6][1])+","+str(stables[num][6][2])+")")

            self.mass3.configure(state="normal")
            self.mass3.delete(0,"end")
            self.mass3.insert(0, stables[num][7])

            self.velocity3.configure(state="normal")
            self.velocity3.delete(0,"end")
            self.velocity3.insert(0, "("+str(stables[num][8][0])+","+str(stables[num][8][1])+","+str(stables[num][8][2])+")")

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

        def update():
            pos = str(self.position1.get())
            if self.position1.get()[0]!="(" or self.position2.get()[0]!="(" or self.position3.get()[0]!="(" or self.position1.get()[-1]!=")" or self.position2.get()[-1]!=")" or self.position3.get()[-1]!=")" or self.velocity1.get()[0]!="(" or self.velocity2.get()[0]!="(" or self.velocity3.get()[0]!="(" or self.velocity1.get()[-1]!=")" or self.velocity2.get()[-1]!=")" or self.velocity3.get()[-1]!=")":
                self.textfield.configure(text="Please include AND close brackets for position and velocity")
                return
            else:
                customData = list((eval(self.position1.get()),eval(self.mass1.get()),eval(self.velocity1.get()),eval(self.position2.get()),eval(self.mass2.get()),eval(self.velocity2.get()),eval(self.position3.get()),eval(self.mass3.get()),eval(self.velocity3.get())))
                return(customData)





        stable1 = list(((0,0,0),10,(0,0,0),(0,0,0),24,(0,0,0),(0,0,0),15,(0,0,0)))
        stable2 = list(((20,0,0),25,(10,0,0),(0,0,10),12,(0,0,20),(10,0,0),30,(20,0,0)))
        stable3 = list(((0,20,0),20,(0,10,0),(0,10,0),6,(0,20,0),(0,10,0),45,(0,20,0)))
        stable4 = list(((0,0,20),20,(0,0,10),(10,0,0),3,(20,0,0),(0,0,10),60,(0,0,20)))
        stables = list((stable1, stable2, stable3, stable4))


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
        self.dropdown1 = customtkinter.CTkOptionMenu(self.dropbox_frame, values=["Stable 1", "Stable 2", "Stable 3", "Stable 4"], command=stableOrbits)
        self.dropdown1.pack(padx=(20,20), pady=20, side="left")

        #button to allow for custom inputs
        var1 = IntVar()
        self.button = customtkinter.CTkCheckBox(master = self.dropbox_frame, text="Custom", variable=var1, onvalue=1, offvalue=0, corner_radius=8, command=custom)
        self.button.pack(padx=20, pady=20, side="left")

        #self.objects = customtkinter.CTkFrame(self)
        #self.objects.grid(column=3, row=0, sticky="nsew", padx=20, pady=20)
        #self.dropdown2 = customtkinter.CTkOptionMenu(self.objects, values=["Object 1", "Object 2", "Object 3"])
        #self.dropdown2.pack(padx=20, pady=20, side="top")


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

        self.submitframe = customtkinter.CTkFrame(self)
        self.submitframe.grid(column=2, row=3, pady=(0,20), padx=20, sticky="nsew")
        self.submit = customtkinter.CTkButton(self.submitframe,text="Simulate", command=update)
        self.submit.pack(pady=20, padx=20,side="top")
        self.textfield = customtkinter.CTkLabel(self.submitframe, height=20, width=300, text="")
        self.textfield.pack(pady=20,side="top")

        #self.variable_frame = customtkinter.CTkFrame(self)
        #self.variable_frame.grid(column=2, row=1, columnspan=2, rowspan=2, padx=20, pady=20, sticky="nsew")
        #self.variable_frame.grid_columnconfigure((0,1), weight=1)
        #self.variable_frame.grid_rowconfigure((0,1), weight=1)
        #self.objects = customtkinter.CTkSegmentedButton(self.variable_frame, values=["Object 1", "Object 2", "Object 3"])
        #self.objects.pack(padx=20, pady=20, side="top")
        #self.positionLabel = customtkinter.CTkLabel(master=self.variable_frame, text="Position (x,y,z):")
        #self.positionLabel.pack(padx=20, pady=(30,5), side="top")
        #self.position = customtkinter.CTkEntry(self.variable_frame, state="disabled")
        #self.position.pack(padx=20, pady=5, side="top")
        #self.massLabel = customtkinter.CTkLabel(master=self.variable_frame, text="Mass (solar mass):")
        #self.massLabel.pack(padx=20, pady=(30,5), side="top")
        #self.mass = customtkinter.CTkEntry(self.variable_frame, state="disabled")
        #self.mass.pack(padx=20, pady=5, side="top")
        #self.velocityLabel = customtkinter.CTkLabel(master=self.variable_frame, text="Velocity (km/s) (x,y,z):")
        #self.velocityLabel.pack(padx=20, pady=(30,5), side="top")
        #self.velocity = customtkinter.CTkEntry(self.variable_frame, state="disabled")
        #self.velocity.pack(padx=20, pady=5, side="top")


        #self.text = customtkinter.CTkLabel(self.variable_frame, text="Velocities (km/s)")
        #self.text.pack(padx=20, pady=20, side="left")
        #self.velocityx = customtkinter.CTkEntry(self.variable_frame, placeholder_text="Velocity (x)")
        #self.velocityx.pack(padx=20, pady=20, side="left")
        #self.velocityy = customtkinter.CTkEntry(self.variable_frame, placeholder_text="Velocity (y)")
        #self.velocityy.pack(padx=20, pady=20, side="left")
        #self.velocityz = customtkinter.CTkEntry(self.variable_frame, placeholder_text="Velocity (z)")
        #self.velocityz.pack(padx=20, pady=20, side="left")





app().mainloop()