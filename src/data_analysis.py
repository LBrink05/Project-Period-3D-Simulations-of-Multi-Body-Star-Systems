
import csv
import numpy as np
from scipy.interpolate import CubicSpline

import scipy

def interpolate_ps(NUM_BODIES, path, precision=None):
    discrete_ps = []
    Interpolated_ps = []

    # Read all data with time index
    for b in range(NUM_BODIES):
        path_b = str(path) + f"/body{b}.csv"
        with open(path_b, "r") as data:
            for i, x in enumerate(data):
                line = x.strip().split(",")
                if len(line) >= 6:
                    discrete_ps.append([
                        b, i,
                        float(line[0]), float(line[1]), float(line[2]),
                        float(line[3]), float(line[4]), float(line[5])
                    ])

    discrete_ps = np.array(discrete_ps)

    # Create spline for each body
    for b in range(NUM_BODIES):
        body_data = discrete_ps[discrete_ps[:, 0] == b]

        if len(body_data) > 1:
            time = body_data[:, 1]              # frame index
            state = body_data[:, 2:]            # x y z vx vy vz

            # Local, non-smoothing interpolation
            spline = CubicSpline(time, state, axis=0)
            Interpolated_ps.append(spline)

    return Interpolated_ps


# Create combined multivariable function of system
def phase_space_function(Interpolated_ps, t):
    result = []
    for spline in Interpolated_ps:
        result.extend(spline(t))  # evaluates [x,y,z,vx,vy,vz]
    return np.array(result)



def difference_function(reference_lines, simulated_lines):
    def diff(t):
        return (
            phase_space_function(reference_lines, t)
            - phase_space_function(simulated_lines, t)
        )
    return diff



def error_function(reference_lines, simulated_lines):
    def percent_error(t):
        ref_ps = phase_space_function(reference_lines, t)
        sim_ps = phase_space_function(simulated_lines, t)

        diff = ref_ps - sim_ps
        ref_norm = np.linalg.norm(ref_ps)

        if ref_norm == 0:
            return float('inf')

        return (np.linalg.norm(diff) / ref_norm) * 100

    return percent_error
