import numpy as np
import time

def create_cylinder_mask(NX, NY, dx, dy, radius, center):
    """Create a mask for a cylinder in a grid.

    Args:
        NX (int): Number of cells in the x direction.
        NY (int): Number of cells in the y direction.
        dx (float): Cell size in the x direction.
        dy (float): Cell size in the y direction.
        radius (float): Radius of the cylinder.
        center (tuple): Center of the cylinder.

    Returns:
        np.ndarray: A mask for the cylinder.
    """
    X, Y = np.meshgrid(np.arange(NX) * dx, np.arange(NY) * dy)
    return (X - center[0])**2 + (Y - center[1])**2 <= radius**2




class Prm:
    def __init__(self):
        self.NX = 100                       # total grid points (including ghost points)
        self.NY = 100                       # total grid points (including ghost points)
        self.nx = self.NX - 2               # real grid points
        self.ny = self.NY - 2               # real grid points
        self.LX = 10.0                      # Length of domain (x direction)
        self.LY = 5.0                       # Length of domain (y direction)
        self.dx = self.LX / self.nx         # grid spacing (x direction)
        self.dy = self.LY / self.ny         # grid spacing (y direction)
        self.dt = 0.001                     # time step
        self.u_init = 1.0                   # entering speed
        self.L = 1.0                        # characteristic length of the object
        self.nu = 1.51e-5                   # kinematic viscosity of air at atmospheric conditions
        self.Re = self.u_init * self.L / self.nu # Reynolds number

        self.U =np.zeros((self.NX, self.NY, 3))       # having u, v and p as components
        
        self.opject_mask = create_cylinder_mask(self.NX, self.NY, self.dx, self.dy, 0.1, (0.5, 0.5))
        
    
    def boundary_conditions(self):
        """
        Boundary conditions at the edge of the domain.
        - on the left side  of the domain, the flow is incoming with a velocity equal to Ue1 and the pressure satisfies ∂p/∂x = 0
        - on the right side of the domain, the flow is free so that ∂u/∂x = ∂v/∂x =0 and p=0,
        - on the horizontal sides, the walls are impenetrable so that v = 0 and a slip condition is imposed so that ∂u/∂y = ∂p/∂y = 0

        Boundary conditions for the object.
        - u = v = 0 and ∂p/∂n  = 0 and slip condition Things are missing here!
        """
        # Applying boundary conditions for the left and right walls
        self.U[0,:,0] = np.full_like(self.U[0,:,0], self.u_init)
        self.U[0,:,1] = np.zeros_like(self.U[0,:,1])
        self.U[0,:,2] = self.U[1,:,2]

        self.U[-1,:,0] = self.U[-2,:,0]
        self.U[-1,:,1] = self.U[-2,:,1]
        self.U[-1,:,2] = np.zeros_like(self.U[-1,:,2])

        # Applying boundary conditions for the top and bottom walls
        self.U[:,0,0] = self.U[:,1,0]  
        self.U[:,0,1] = np.zeros_like(self.U[:,0,1])  
        self.U[:,0,2] = self.U[:,1,2]  

        self.U[:,-1,0] = self.U[:,-2,0]
        self.U[:,-1,1] = np.zeros_like(self.U[:,-1,1])
        self.U[:,-1,2] = self.U[:,-2,2]

        # Boundary conditions for the object
        self.U[self.opject_mask, 0] = 0
        self.U[self.opject_mask, 1] = 0
        '''HERE THE PRESSURE GRADIENT IS MISSING'''
        self.U[self.opject_mask, 2] = 0
        self.U[self.opject_mask, 2] = 0
    
        
    def Semilag(self, U):
        sign_u = np.sign(U[:,:,0])
        sign_v = np.sign(U[:,:,1])
        aux = np.zeros((self.NX, self.NY))

        a = 1 - sign_u * self.dt / self.dx
        b = 1 - sign_v * self.dt / self.dy
        a_positive = np.where(sign_u > 0, a, 1 / a)
        b_positive = np.where(sign_v > 0, b, 1 / b)

        q = U[:,:,2]
        q0 = q[1:-1, 1:-1]
        q1 = q[np.maximum(sign_u[1:-1, 1:-1], 0), np.maximum(sign_v[1:-1, 1:-1], 0)]
        q2 = q[np.maximum(sign_u[1:-1, 1:-1], 0), np.minimum(sign_v[1:-1, 1:-1], self.NY - 1)]
        q3 = q[np.minimum(sign_u[1:-1, 1:-1], self.NX - 1), np.maximum(sign_v[1:-1, 1:-1], 0)]

        aux = (a_positive * b_positive * q0 +
            (1 - a_positive) * b_positive * q1 +
            a_positive * (1 - b_positive) * q2 +
            (1 - a_positive) * (1 - b_positive) * q3)

        return aux



"""
def Semilag(u, v, q, prm, sign):
    sign_u = np.sign(sign * u)
    sign_v = np.sign(sign * v)
    aux = np.zeros((prm.NX, prm.NY))
    for i in range(1, prm.NX - 1):
        for j in range(1, prm.NY - 1):
            a = 1 - sign_u[i, j] * prm.dt / prm.dx if sign_u[i, j] > 0 else 1 + sign_u[i, j] * prm.dt / prm.dx
            b = 1 - sign_v[i, j] * prm.dt / prm.dy if sign_v[i, j] > 0 else 1 + sign_v[i, j] * prm.dt / prm.dy
            aux[i, j] = (a * b * q[i, j] +
                         (1 - a) * b * q[i - sign_u[i, j], j] +
                         a * (1 - b) * q[i, j - sign_v[i, j]] +
                         (1 - a) * (1 - b) * q[i - sign_u[i, j], j - sign_v[i, j]])
    return aux

def Semilag2(u, v, q, prm):
    q0 = q.copy()
    q = Semilag(u, v, q, prm, 1)
    q = Semilag(u, v, q, prm, -1)
    q = q0 + (q0 - q) / 2
    return Semilag(u, v, q, prm, 1)

def boundary_conditions(u, v, p, prm):
    u[0, :] = 2 * prm.U - u[1, :]
    u[-1, :] = u[-2, :]
    v[0, :] = -v[1, :]
    v[-1, :] = -v[-2, :]
    p[0, :] = p[1, :]
    p[-1, :] = -p[-2, :]
    u[:, 0] = u[:, 1]
    u[:, -1] = u[:, -2]
    v[:, 0] = -v[:, 1]
    v[:, -1] = -v[:, -2]
    p[:, 0] = p[:, 1]
    p[:, -1] = p[:, -2]

def init(prm):
    prm.NX = 100                         # total grid points (including ghost points)
    prm.NY = 100                         # total grid points (including ghost points)
    prm.nx = prm.NX - 2                 # real grid points
    prm.ny = prm.NY - 2                 # real grid points
    prm.LX = 10.0                        # Length of domain (x direction)
    prm.LY = 5.0                        # Length of domain (y direction)
    prm.dx = prm.LX / prm.nx          # grid spacing (x direction)
    prm.dy = prm.LY / prm.ny          # grid spacing (y direction)
    prm.dt = 0.001                      # time step
    prm.U = 1.0                          # entering speed
    prm.L = 1.0                          # characteristic length of the object
    prm.nu = 1.51e-5                    # kinematic viscosity of air at atmospheric conditions
    prm.Re = prm.U * prm.L / prm.nu   # Reynolds number

prm = Prm()
init(prm)

t = 0.0

# allocate memory
u = np.zeros((prm.NX, prm.NY))  # x component of velocity, initialized to 0
v = np.zeros((prm.NX, prm.NY))  # y component of velocity, initialized to 0
p = np.zeros((prm.NX, prm.NY))  # pressure, initialized to 0

# count time
start = time.time()
while t < 2:
    # advection term
    u = Semilag2(u, v, u, prm)
    v = Semilag2(u, v, v, prm)
    t += prm.dt
end = time.time()
print("time =", end - start)

"""