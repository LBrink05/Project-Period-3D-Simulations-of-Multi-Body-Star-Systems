#lagrange triangle solution
NAMES = ["Lagrange-triangle-solution", "Triangle-solution", "Triangle", "0"]
data = [["NUM_BODIES", "GRAVITY", "COLORS", "MASS", "START POSITION", "START VELOCITY"]]

#default settings
#CONSTANTS
NUM_BODIES = 3
GRAV = 1
COLOR = [RED,GREEN,BLUE]
MASS = np.array([1,1,1])

#Solution specific configurations

r = 5.0  # distance from center of mass to each body
R = r * np.sqrt(3)
v_mag = np.sqrt(GRAV * MASS[0] / R)

# 120-degree separation
START_POS = np.array([[ r, 0, 0],[-0.5*r,  np.sqrt(3)/2*r, 0],[-0.5*r, -np.sqrt(3)/2*r, 0]], dtype=float)

START_VEL = np.array([[0, v_mag, 0],[-v_mag*np.sqrt(3)/2, -v_mag/2, 0],[ v_mag*np.sqrt(3)/2, -v_mag/2, 0]], dtype=float)


