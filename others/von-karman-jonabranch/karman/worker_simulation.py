import numpy as np
import h5py
import copy
from tqdm import tqdm
from numba import jit, prange
import time
# own package defining the initial conditions and boundary conditions
from helper_setup import *

# related to the fluxes
@jit(nopython=True)
def pressure(rho,rhou,rhov,rhoE, gamma):
    """
    Calculates the pressure of the fluid using the state variables.

    Args:
        rho (float): Density of the fluid.
        rhou (float): Momentum density in the x-direction.
        rhov (float): Momentum density in the y-direction.
        rhoE (float): Total energy density.

    Returns:
        float: Pressure of the fluid.
    """
    
    return (gamma-1)*(rhoE-0.5*(rhou**2+rhov**2)/rho)

@jit(nopython=True)
def char_speed(p, rho, gamma):
    """
    Calculates the characteristic speed of the fluid.

    Args:
        p (float): Pressure of the fluid.
        rho (float): Density of the fluid.

    Returns:
        float: Characteristic speed of the fluid.
    """
    
    return np.sqrt(gamma*p/rho)


@jit(nopython=True)
def get_adaptive_dt(U, NX, NY, dx, gamma, cfl):
    """
    Calculates the time step for the simulation.

    Args:
        U (ndarray): The state of the system.
        dx (float): The spatial step size.
        gamma (float): The specific heat ratio.
        cfl (float): The Courant-Friedrichs-Lewy number.

    Returns:
        float: The time step.
    """
    # Unpack the U array into separate variables for better performance
    rho = U[:,:,0]
    u = U[:,:,1] / rho
    v = U[:,:,2] / rho
    p = (gamma - 1) * (U[:,:,3] - 0.5 * rho * (u**2 + v**2))
    c = np.sqrt(gamma * p / rho)

    # Compute the local time step using vectorized operations
    dt_ij = cfl * dx / (np.sqrt(u**2 + v**2) + c)

    # Find the minimum value of dt_ij
    min_dt = np.min(dt_ij)

    # Check for negative or zero pressure
    if np.any(p <= 0):
        raise ValueError("Negative or zero pressure detected.")

    return min_dt

@jit(nopython=True)
def fluxFunction(Ul, Ur,mark_r, mark_l, direction, gamma):
    """
    Calculates the flux between two adjacent cells.

    Args:
        Ul (ndarray): State vector of the left cell.
        Ur (ndarray): State vector of the right cell.
        direction (int): Direction of the flux. 0 for x-direction, 1 for y-direction.

    Returns:
        ndarray: Flux between the Ul and Ur cells.
    """

    if mark_r and direction == 0:
        Ur[0] = Ul[0]
        Ur[3] = Ul[3]
        Ur[1] = -Ul[1]
        Ur[2] = 0
    elif mark_r and direction == 1:
        Ur[0] = Ul[0]
        Ur[3] = Ul[3]
        Ur[1] = 0
        Ur[2] = -Ul[2]
    elif mark_l and direction == 0:
        Ul[0] = Ur[0]
        Ul[3] = Ur[3]
        Ul[1] = -Ur[1]
        Ul[2] = 0        
    elif mark_r and direction == 1:
        Ul[0] = Ur[0]
        Ul[3] = Ur[3]
        Ul[1] = 0
        Ul[2] = -Ur[2]
    # pressure and sound speed
    pl = pressure(Ul[0],Ul[1],Ul[2], Ul[3], gamma)
    pr = pressure(Ur[0],Ur[1],Ur[2], Ur[3], gamma)
    cl = char_speed(pl, Ul[0], gamma)
    cr = char_speed(pr, Ur[0], gamma)        

    # characteristic speed
    if Ul[0]*cl > Ur[0]*cr:
        a = 1.1*Ul[0]*cl
    else:
        a = 1.1*Ur[0]*cr
    

    if direction == 0:
        ustar = 0.5*(Ul[1]/Ul[0] + Ur[1]/Ur[0]) - 0.5*(pr - pl)/a
        pstar = 0.5*(pl + pr) - 0.5*(Ur[1]/Ur[0]-Ul[1]/Ul[0])*a

        if ustar >= 0:
            Ubar = Ul
        else:
            Ubar = Ur
    
        return np.array([Ubar[0]*ustar, Ubar[1]*ustar + pstar, Ubar[2]*ustar, Ubar[3]*ustar + pstar*ustar])
    else :
        ustar = 0.5*(Ul[2]/Ul[0] + Ur[2]/Ur[0]) - 0.5*(pr - pl)/a
        pstar = 0.5*(pl + pr) - 0.5*(Ur[2]/Ur[0]-Ul[2]/Ul[0])*a

        if ustar >= 0:
            Ubar = Ul
        else:
            Ubar = Ur
        return np.array([Ubar[0]*ustar, Ubar[1]*ustar, Ubar[2]*ustar + pstar, Ubar[3]*ustar + pstar*ustar])
    


@jit(nopython=True, parallel=True)
def update_flux(Uold, U, mark, dt, NX, NY, dx, gamma):
    # Update the solution using the indices computed above
    for i in prange(NX):  # Parallelize the outer loop
        for j in range(NY):
            # Compute the indices for the neighboring cells with periodic boundary conditions
            i_plus_one = (i + 1) % NX
            i_minus_one = (i - 1) % NX
            j_plus_one = (j + 1) % NY
            j_minus_one = (j - 1) % NY
            
            # Update the solution using the indices computed above
            U[i, j] = Uold[i, j] - dt/dx * (
                fluxFunction(Uold[i, j], Uold[i_plus_one, j], mark[i, j], mark[i_plus_one, j], 0, gamma) 
                - fluxFunction(Uold[i_minus_one, j], Uold[i, j], mark[i_minus_one, j], mark[i, j], 0, gamma) 
                + fluxFunction(Uold[i, j], Uold[i, j_plus_one], mark[i, j], mark[i, j_plus_one], 1, gamma) 
                - fluxFunction(Uold[i, j_minus_one], Uold[i, j], mark[i, j_minus_one], mark[i, j], 1, gamma)
            )
            #Finding errors
            if U[i,j,0] <= 0 or U[i,j,3] <= 0:
                print(f"Error low at {i, j}\n")
            if U[i,j,0] > 10000 or U[i,j,3] > 10000:
                print(f"Error high at {i, j}\n")
    return U
# Simulation

def run_simulation(params):
    """
    Run the simulation.

    Args:
        params: Dictionary containing the simulation parameters.
    """
    
    # Extract parameters
    initial_condition_name = params["initial_conditions"]
    initial_condition = globals()[initial_condition_name]
    boundary_conditions_name = params["boundary_conditions"]
    boundary_conditions = globals()[boundary_conditions_name]
    object_condition_name = params["object_condition"]
    object_condition = globals()[object_condition_name]
    t_start = params["t_start"]
    t_end = params["t_end"]
    dt_given = params["dt"]
    dt_io = params["dt_io"]
    filename_save = params["filename_save"]
    dx = params["dx"]
    xmin = params["xmin"]
    xmax = params["xmax"]
    ymin = params["ymin"]
    ymax = params["ymax"]
    gamma = params["gamma"]
    cfl = params["cfl"]

    
    # Initialize the 2D grid
    X, Y = np.arange(xmin, xmax, dx), np.arange(ymin, ymax, dx)
    NX, NY = len(X), len(Y)
    U_init = initial_condition(X, Y)
    mark = object_condition(X,Y)
    U = copy.deepcopy(U_init)
    ts = [t_start]
    start_time = time.time()
    # Run the simulation
    with h5py.File(filename_save, "w") as f:
        f.create_dataset("X", data=X)
        f.create_dataset("Y", data=Y)
        f.create_dataset(f"U_{ts[-1]}", data=U)

        # Run simulation
        t_last = t_start
        progress_bar = tqdm(total=(t_end - t_start))
        while ts[-1] < t_end:
            Uold = copy.deepcopy(U)

            Uold = boundary_conditions(Uold)

            if dt_given == 0:
                dt = get_adaptive_dt(Uold, NX, NY, dx, gamma, cfl)
            else:
                dt = dt_given

            if ts[-1] + dt >= t_end:
                dt = t_end - ts[-1]
            
            ts.append(ts[-1] + dt)

            progress_bar.update(round(dt, 5))

            # Update the solution using the indices computed above
            U = update_flux(Uold, U, mark, dt, NX, NY, dx, gamma)

            if dt_io == 0:
                f.create_dataset(f"U_{ts[-1]}", data=U)
                t_last = ts[-1]
            elif ts[-1] - t_last >= dt_io or ts[-1] == t_end:
                f.create_dataset(f"U_{ts[-1]}", data=U)
                t_last = ts[-1]
            
        progress_bar.close()
    print(f"Solution saved to: {filename_save}")
    end_time = time.time()
    return end_time - start_time

