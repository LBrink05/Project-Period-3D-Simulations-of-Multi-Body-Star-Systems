import numpy as np
import csv

#Used to calculate the stable conditions using spherical coordinates.

# Convert to Python floats immediately
rt32 = float(np.sqrt(3)/2)
v = float(np.sqrt(1/(5*np.sqrt(3))))
rt32_times_5 = float(rt32 * 5)
v_times_rt32 = float(v * rt32)
v_div_2 = float(v / 2)

# Create the list with pre-computed float values
stable1 = [
    (5, 0, 0),
    1.0,
    (0.0, v, 0.0),
    (-2.5, rt32_times_5, 0.0),
    1.0,
    (-v_times_rt32, -v_div_2, 0.0),
    (-2.5, -rt32_times_5, 0.0),
    1.0,
    (v_times_rt32, -v_div_2, 0.0),
    "Equilateral Triangle"
]

print(stable1)


'''stable1 = list(((5, 0, 0),1,(0,v,0),(-0.5*5,rt32*5,0),1,(-v*rt32,-v/2,0),(-0.5*5,-rt32*5,0),1,(v*rt32,-v/2,0)))
stable2 = list(((2.57429,0,0),1,(0.216343, 0.332029,0),(-2.57429,0,0),1,(0.216343, 0.332029,0),(0,0,0),1,(-0.432686, -0.664058,0)))
stable3 = list(((1,0,0),1,(0,0.5, 0),(-0.5,rt32,0),1,(-0.433,-0.25, 0),(-0.5,-rt32,0),1,(0.433,-0.25, 0)))
stable4 = list((( 1.0, 0.0, 0), 1, ( 0.0,  omega, 0),(-0.5,  rt32, 0), 1, (-omega*rt32, -omega/2, 0),(-0.5, -rt32, 0), 1, ( omega*rt32, -omega/2, 0)))
stable5 = list(((-13.02117, 0, 0), 1, ( 0.085236,  0.03498, 0),( 13.02117, 0, 0), 1, ( 0.085236,  0.03498, 0),( 0,        0, 0), 1, (-0.170472, -0.06996, 0))) #butterfly I
stable6 = list(((-12.04917, 0, 0), 1, ( 0.113643,  0.028219, 0),( 12.04917, 0, 0), 1, ( 0.113643,  0.028219, 0),( 0,        0, 0), 1, (-0.227286, -0.056438, 0)))  #butteryfly II
stable7 = list(((-8.29368, 0, 0), 1, ( 0.161277,  0.137507, 0),( 8.29368, 0, 0), 1, ( 0.161277,  0.137507, 0),( 0,       0, 0), 1, (-0.322554, -0.275014, 0))) #moth I
stable8 = list(((-7.83516, 0, 0), 1, ( 0.156829,  0.16207, 0),( 7.83516, 0, 0), 1, ( 0.156829,  0.16207, 0),( 0,       0, 0), 1, (-0.313658, -0.32414, 0))) #moth II
stable9 = list(((-7.17921, 0, 0), 1, ( 0.208677,  0.130401, 0),( 7.17921, 0, 0), 1, ( 0.208677,  0.130401, 0),( 0,       0, 0), 1, (-0.417354, -0.260802, 0))) #Yarn
stable10 = list(((-8.57406, 0, 0), 1, ( 0.175521,  0.104039, 0),( 8.57406, 0, 0), 1, ( 0.175521,  0.104039, 0),( 0,       0, 0), 1, (-0.351042, -0.208078, 0))) #yinyang
stables = list((stable1, stable2, stable3, stable4, stable5, stable6, stable7, stable8, stable9, stable10))'''