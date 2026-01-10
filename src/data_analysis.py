import numpy as np
import csv
import sklearn
import scipy

def interpolate_data(NUM_BODIES, path, precision):

    trails = []
    Interpolated_trails = []

    # Read all data with time index
    for b in range(NUM_BODIES):
        path_b = str(path) + f"/body{b}.csv"
        with open(path_b, "r") as data:
            for i, x in enumerate(data):
                line = x.strip().split(",")
                if line and len(line) >= 3:
                    trails.append([b, i, float(line[0]), float(line[1]), float(line[2])])

    trails = np.array(trails)

    # Create continuous function for each body
    for b in range(NUM_BODIES):
        body_data = trails[trails[:, 0] == b]
        
        if len(body_data) > 1:
            # X = time, y = [x, y, z] positions
            time = body_data[:, 1].reshape(-1, 1)
            positions = body_data[:, 2:]  # x, y, z
            
            # Fit the model
            model = LinearRegression()
            model.fit(time, positions)
            
            # Store the model itself - this is your continuous function
            Interpolated_trails.append(model)

    return Interpolated_trails #continous time function

# Now you can evaluate at ANY time
def get_continous_position(body_index, time):
    #Get position of body at any time
    return Interpolated_trails[body_index].predict([[time]])[0]


def difference_function(body_index, reference_lines, simulated_lines):
    #Returns a function that gives the difference between two bodies' positions
    func1 = reference_lines[body_index]
    func2 = simulated_lines[body_index]
    
    def diff(t):
        return func1(t) - func2(t)
    
    return diff

def integrate_trajectory(interpolated_trail ,body_index, t_start, t_end):
    #Integrate position over time for a body
    interp_func = Interpolated_trails[body_index]
    
    # Integrate each coordinate separately
    x_integral, _ = quad(lambda t: interp_func(t)[0], t_start, t_end)
    y_integral, _ = quad(lambda t: interp_func(t)[1], t_start, t_end)
    z_integral, _ = quad(lambda t: interp_func(t)[2], t_start, t_end)
    
    return np.array([x_integral, y_integral, z_integral])

    # Example: integrate body 0 from time 0 to 10
    result = integrate_trajectory(0, 0, 10)
    print(f"Integral: x={result[0]}, y={result[1]}, z={result[2]}")

def percent_error(reference_lines, simulated_lines):


     return percent_error