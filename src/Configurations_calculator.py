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
