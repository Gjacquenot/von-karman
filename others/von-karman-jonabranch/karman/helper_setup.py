import numpy as np

def prim2cons(rho, u, v, p, gamma):
    """
    Converts the primitive variables to the conservative variables.

    Args:
        rho (float): Density of the fluid.
        u (float): Velocity in the x-direction.
        v (float): Velocity in the y-direction.
        p (float): Pressure of the fluid.

    Returns:
        ndarray: State vector in conservative variables.
    """
    return rho, rho*u, rho*v, p/(gamma-1) + 0.5*rho*(u**2+v**2)

def cons2prim(rho, rhou, rhov, rhoE, gamma):
    """
    Converts the conservative variables to the primitive variables.

    Args:
        U (ndarray): State vector in conservative variables.

    Returns:
        float: Density of the fluid.
        float: Velocity in the x-direction.
        float: Velocity in the y-direction.
        float: Pressure of the fluid.
    """
    u = rhou/rho
    v = rhov/rho
    p = (gamma-1)*(rhoE - 0.5*(rhou**2+rhov**2)/rho)
    return rho, u, v, p




#-----------------# Boundary conditions #-----------------#
def kelvin_bc(U):
    """
    Apply periodic boundary conditions left right, neumann top bottom

    Args:
        U (ndarray): State vector.

    Returns:
        ndarray: State vector with boundary conditions applied.
    """
    U[:, 0, :] = U[:, 1, :]
    U[:, -1, :] = U[:, -2, :]
    return U

def periodic_bc(U):
    """
    Apply periodic boundary conditions

    Args:
        U (ndarray): State vector.

    Returns:
        ndarray: State vector with boundary conditions applied.
    """
    return U

def neumann_bc(U):
    """
    Apply periodic boundary conditions

    Args:
        U (ndarray): State vector.

    Returns:
        ndarray: State vector with boundary conditions applied.
    """
    U[0, :, :] = U[1, :, :]
    U[-1, :, :] = U[-2, :, :]
    U[:, 0, :] = U[:, 1, :]
    U[:, -1, :] = U[:, -2, :]
    return U

def karman_bc(U):
    """
    Apply boundary condition for karman vortex street

    Args:
        U (ndarray): State vector.

    Returns:
        ndarray: State vector with boundary conditions applied.
    """   

    U[0,:,0], U[0,:,1], U[0,:,2], U[0,:,3] = U[1,:,0], 0.1, 0, U[0,:,3]
    U[-1,:,0], U[-1,:,1], U[-1,:,2], U[-1,:,3] = U[-2,:,0], U[-2,:,1], U[-2,:,2], U[-2,:,3]
    U[:,0,0], U[:,0,1], U[:,0,2], U[:,0,3] = U[:,1,0], 0, -U[:,1, 2] , U[:,1,3]
    U[:,-1,0], U[:,-1,1], U[:,-1,2], U[:,-1,3] = U[:,-2,0], 0, -U[:,-2, 2] , U[:,-2,3]
    return U
    

#-----------------# Initial conditions #-----------------#
def initial_conditions_turned(X, Y):
    """
    Generate initial conditions for a problem where the conditions change based on the Y coordinate.

    Args:
        X (ndarray): Array of X coordinates.
        Y (ndarray): Array of Y coordinates.

    Returns:
        ndarray: Initial state vector for the problem with Y-based initial conditions.
    """
    NX, NY = len(X), len(Y)    
    U_init = np.zeros((NX, NY, 4))
    for i in range(NX):
        for j in range(NY):
            x, y = X[i], Y[j]
            if x+y < 1:
                U_init[i, j, 0] = 1  
                U_init[i, j, 3] = 1
            else:
                U_init[i, j, 0] = 0.125  
                U_init[i, j, 3] = 0.1
    return U_init

#DEAL WITH THIS PRIM TO CONS
def initial_kelvin_helmholtz(X, Y):
    """
    Generate initial conditions for the Kelvin-Helmholtz instability problem using vectorized operations.

    Args:
        X (ndarray): Array of X coordinates.
        Y (ndarray): Array of Y coordinates.

    Returns:
        ndarray: Initial state vector for the Kelvin-Helmholtz instability problem.
    """
    y1 = 0.5
    y2 = 1.5
    a = 0.05
    amp = 0.01
    sigma = 0.2

    X_grid, Y_grid = np.meshgrid(X, Y, indexing='ij')
    
    U_init = np.zeros((X_grid.shape[0], Y_grid.shape[1], 4))

    y_term = (Y_grid - y1) / a
    rho, rhou, rhov, rhoE = prim2cons(1 + np.tanh(y_term) - np.tanh((Y_grid - y2) / a),
                                    np.tanh(y_term) - np.tanh((Y_grid - y2) / a),
                                    amp * np.sin(np.pi * X_grid) * (np.exp(- (Y_grid - y1) ** 2 / (sigma ** 2)) 
                                                              + np.exp(- (Y_grid - y2) ** 2 / (sigma ** 2))),
                                    10, 1.4)                               

    U_init[:,:,0] = rho
    U_init[:,:,1] = rhou
    U_init[:,:,2] = rhov
    U_init[:,:,3] = rhoE
    
    return U_init


def initial_2dR(X, Y):#LEADS TO PROBLEMS IN THE SOLVER
    """
    Generate initial conditions for the blast problem using vectorized operations.

    Args:
        X (ndarray): Array of X coordinates.
        Y (ndarray): Array of Y coordinates.

    Returns:
        ndarray: Initial state vector for the Kelvin-Helmholtz instability problem.
    """
    X_grid, Y_grid = np.meshgrid(X, Y, indexing='ij')
    
    U_init = np.full((X_grid.shape[0], Y_grid.shape[1], 4),0.0001)
    mask1 = (X_grid > 0.8) & (Y_grid > 0.8)
    mask2 = (X_grid <= 0.8) & (Y_grid > 0.8)
    mask3 = (X_grid <= 0.8) & (Y_grid <= 0.8)
    mask4 = (X_grid > 0.8) & (Y_grid <= 0.8)
    U_init[mask1] = [prim2cons(1.5, 0., 0., 1.5, 1.4)]
    U_init[mask2] = [prim2cons(0.5323, 1.206, 0., 0.3, 1.4)]
    U_init[mask3] = [prim2cons(0.1380, 1.206, 1.206, 0.029, 1.4)]
    U_init[mask4] = [prim2cons(0.5323, 0., 1.206, 0.3, 1.4)]
    
    return U_init

def initial_sod(X, Y):
    """
    Generate initial conditions for the blast problem using vectorized operations.

    Args:
        X (ndarray): Array of X coordinates.
        Y (ndarray): Array of Y coordinates.

    Returns:
        ndarray: Initial state vector for the Kelvin-Helmholtz instability problem.
    """
    X_grid, Y_grid = np.meshgrid(X, Y, indexing='ij')
    
    U_init = np.zeros((X_grid.shape[0], Y_grid.shape[1], 4))
    mask1 = (X_grid <= 0.5)
    mask2 = (X_grid >= 0.5) 

    U_init[mask1] = [1., 0., 0., 1.]
    U_init[mask2] = [0.125, 0., 0., 0.1]
 

    return U_init

def initial_karman(X, Y):
    """
    Generate initial conditions for the karman vortex street
    
    Args:
        X (ndarray): Array of X coordinates.
        Y (ndarray): Array of Y coordinates.

    Returns:
        ndarray: Initial state vector for the Kelvin-Helmholtz instability problem.
    """
    X_grid, Y_grid = np.meshgrid(X, Y, indexing='ij')
    
    U_init = np.zeros((X_grid.shape[0], Y_grid.shape[1], 4))
    U_init[1,:,0] = 1
    U_init[:,:,0] = 20
    U_init[:,:,3] = 20

    return U_init

#---------------------# OBJECTS AND OTHER STUFF #---------------------#

def square(X, Y):
    X_grid, Y_grid = np.meshgrid(X, Y, indexing='ij')
    mask = (X_grid > 0.4) & (X_grid < 0.6) & (Y_grid > 0.4) & (Y_grid < 0.6)
    return mask

def noobject(X, Y):
    X_grid, Y_grid = np.meshgrid(X, Y, indexing='ij')
    mask = (X_grid > 1) & (X_grid < 0.6)
    return mask