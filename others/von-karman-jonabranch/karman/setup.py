# Simulation Setup

# Function for initial conditions
# Choose one of: initial_kelvin_helmholtz, initial_2dR, initial_sod, initial_conditions_turnedm initial_karman
initial_conditions = "initial_karman"

# Function for boundary conditions
# Choose one of: periodic_bc, neumann_bc, kelvin_bc, karman_bc
boundary_conditions = "karman_bc"

# Function to describe if there is an object in the domain
# Choose one of: noobject, square
object_condition = "square"

# Time parameters
t_start = 0      # Start time
t_end = 100       # End time
dt = 0.0           # Time step (if zero, the CFL condition is used)
dt_io = 1     # Time step for input/output

# File and naming
name = "karman_t100"                    # Name for the simulation
filename_save = f"results/{name}.h5"            # File to save simulation results

# Spatial parameters
dx = 0.01        # Grid step
xmin, xmax = 0, 3    # Minimum and maximum X coordinates
ymin, ymax = 0, 1    # Minimum and maximum Y coordinates

# Fluid properties
gamma = 1.4      # Ratio of specific heats (for compressible flow)

# CFL condition
cfl = 0.4        # CFL number for stability (only used if dt = 0)
