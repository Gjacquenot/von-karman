from src import pressure_poisson, derive
import numpy as np
import h5py
import copy
from tqdm import tqdm
from numba import jit, prange
import time
import os
from src import ip_op
# own package defining the initial conditions and boundary conditions
# from src.helper_setup import *

from datetime import datetime


def BC(u, v, p):

    u[:, 0] = 0
    u[:, -1] = 0
    u[0, :] = 1
    u[-1, :] = u[-2, :]

    v[0, :] = 0
    v[-1, :] = v[-2, :]
    v[:, 0] = 0
    v[:, -1] = 0

    p[0, :] = p[1, :]
    p[-1, :] = p[-2, :]
    p[:, 0] = p[:, 1]
    p[:, -1] = p[:, -2]

    return u, v, p


# -----------------# adaptive timestep #-----------------#
def compute_dt(U, dx, CFL):
    """
    Compute the time step size based on the CFL condition.

    Args:
        U: The solution at the current time step.
        dx: Grid spacing.
        CFL: The Courant-Friedrichs-Lewy number.

    Returns:
        dt: The time step size.
    """
    # Extract components of the solution
    u = U[:, :, 0]
    v = U[:, :, 1]
    dt = CFL * dx**2 / (np.max(np.abs(u)) + np.max(np.abs(v)))
    return dt
# -----------------# Advection #-----------------#


def Semilag(u, v, q, dx, dy, dt, NX, NY):
    """
    1st order semi-Lagrangian advection
    """
    ADVq = np.zeros_like(q)

    # Central domain indices (excluding boundaries)
    i_min, i_max = 1, NX - 1  # Exclude the first and last row
    j_min, j_max = 1, NY - 1  # Exclude the first and last column

    # Matrices where 1 is right, 0 is left or center
    Mx2 = np.sign(np.sign(u[i_min:i_max, j_min:j_max]) + 1.)
    Mx1 = 1. - Mx2

    # Matrices where 1 is up, 0 is down or center
    My2 = np.sign(np.sign(v[i_min:i_max, j_min:j_max]) + 1.)
    My1 = 1. - My2

    # Matrices of absolute values for u and v
    au = abs(u[i_min:i_max, j_min:j_max])
    av = abs(v[i_min:i_max, j_min:j_max])

    # Matrices of coefficients respectively central, external, same x, same y
    Cc = (dx - au * dt) * (dy - av * dt) / dx / dy
    Ce = dt * dt * au * av / dx / dy
    Cmx = (dx - au * dt) * av * dt / dx / dy
    Cmy = dt * au * (dy - dt * av) / dx / dy

    # Computes the advected quantity for the interior (central domain)
    ADVq[i_min:i_max, j_min:j_max] = (Cc * q[i_min:i_max, j_min:j_max] +
                                      Ce * (Mx1 * My1 * q[i_min + 1:i_max + 1, j_min + 1:j_max + 1] +
                                            Mx2 * My1 * q[i_min - 1:i_max - 1, j_min + 1:j_max + 1] +
                                            Mx1 * My2 * q[i_min + 1:i_max + 1, j_min - 1:j_max - 1] +
                                            Mx2 * My2 * q[i_min - 1:i_max - 1, j_min - 1:j_max - 1]) +
                                      Cmx * (My1 * q[i_min:i_max, j_min + 1:j_max + 1] +
                                             My2 * q[i_min:i_max, j_min - 1:j_max - 1]) +
                                      Cmy * (Mx1 * q[i_min + 1:i_max + 1, j_min:j_max] +
                                             Mx2 * q[i_min - 1:i_max - 1, j_min:j_max]))

    q[1:-1, 1:-1] = ADVq[1:-1, 1:-1]  # to keep the boundary conditions
    return q


####
def Semilag2(u, v, q, dx, dy, dt, NX, NY):
    """
    Second order semi-Lagrangian advection
    """

    ADVq = Semilag(u, v, q, dx, dy, dt, NX, NY)
    aux = Semilag(-u, -v, ADVq, dx, dy, dt, NX, NY)
    ADVq = q + (q - aux) / 2.
    ADVq = Semilag(u, v, ADVq, dx, dy, dt, NX, NY)

    return ADVq


@jit(nopython=True, parallel=True)
def upwind_advection(U, V, phi, dt, dx, dy, NX, NY):
    """
    Perform advection using an explicit upwind scheme.

    Args:
        U: The velocity field in the x-direction.
        V: The velocity field in the y-direction.
        phi: The scalar quantity to be advected.
        dt: Time step size.
        dx, dy: Grid spacing in the x and y directions.
        NX, NY: Number of grid points in the x and y directions.

    Returns:
        phi_new: The advected scalar field.
    """

    phi_new = np.zeros_like(phi)
    for i in prange(1, NX - 1):
        for j in range(1, NY - 1):
            u = U[i, j]
            v = V[i, j]

            # Upwind differencing for advection in x-direction
            if u >= 0:
                phi_x = (phi[i, j] - phi[i - 1, j]) / dx
            else:
                phi_x = (phi[i + 1, j] - phi[i, j]) / dx

            # Upwind differencing for advection in y-direction
            if v >= 0:
                phi_y = (phi[i, j] - phi[i, j - 1]) / dy
            else:
                phi_y = (phi[i, j + 1] - phi[i, j]) / dy

            # Update phi using upwind differencing
            phi_new[i, j] = phi[i, j] - dt * (u * phi_x + v * phi_y)

    # Apply boundary conditions for phi_new as needed
    phi[1:-1, 1:-1] = phi_new[1:-1, 1:-1]
    # print(phi[0,0])
    return phi

# -----------------# Diffusion #-----------------#


def compute_diffusion(u, dx, dy):
    # Initialize Laplacian arrays
    lap_u = np.zeros_like(u)

    # Laplacian of u
    lap_u[:, 1:-1] = (u[:, 2:] - 2 * u[:, 1:-1] + u[:, :-2]) / dx**2
    lap_u[1:-1, :] += (u[2:, :] - 2 * u[1:-1, :] + u[:-2, :]) / dy**2

    return lap_u


# -----------------# pressure #-----------------#


def update_pressure(ustarstar, vstarstar, rho, epsilon,
                    dx, dy, nx_sp, ny_sp, K, dt, bc):
    # next compute the pressure RHS: prhs = div(un)/dt + div( [urhs, vrhs])
    prhs = rho * derive.div((1 - epsilon) * ustarstar,
                            (1 - epsilon) * vstarstar, dx, dy, bc) / dt

    # p, err = pressure_poisson.solve_new(p, dx, dy, prhs)

    p = pressure_poisson.solve_spectral(nx_sp, ny_sp, K, prhs)

    return p


def update_velocities_with_pressure(ustarstar, vstarstar, p, dx, dy, dt, bc):
    # Step4
    # finally compute the true velocities
    # u_{n+1} = uh - dt*dpdx
    dpdx = derive.ddx(p, dx, bc)

    dpdy = derive.ddy(p, dy, bc)

    u = ustarstar - dt * dpdx
    v = vstarstar - dt * dpdy

    return u, v, dpdx, dpdy

# -----------------# putting all together for a step #-----------------#


def update_step(u, v, p, dt, NX, NY, dx, epsilon, K, bc):
    """
    Update the solution using the chorin method and semi-Lagrangian method.

    Args:
        Uold: The solution at the previous time step.
        U: The solution at the current time step.
        dt: The time step size.
        NX: Number of grid points in the x-direction.
        NY: Number of grid points in the y-direction.
        dx: Grid spacing.

    Returns:
        U: The updated solution.
    """
    rho = 1

    u, v, p = BC(u, v, p)
    # Advect the solution using the semi-Lagrangian method
    u = Semilag2(u, v, u, dx, dx, dt, NX, NY)
    v = Semilag2(u, v, v, dx, dx, dt, NX, NY)
    u, v, p = BC(u, v, p)
    # Compute diffusion term
    lap_u = compute_diffusion(u, dx, dx)
    lap_v = compute_diffusion(v, dx, dx)
    u += dt * lap_u
    v += dt * lap_v
    u, v, p = BC(u, v, p)
    u, v, p, epsilon = u.T, v.T, p.T, epsilon.T
    # Compute pressure term
    p = update_pressure(u, v, rho, epsilon, dx, dx, NX, NY, K, dt, bc)

    # Update velocity field
    u, v, dpdx, dpdy = update_velocities_with_pressure(u, v, p, dx, dx, dt, bc)
    u, v, p, epsilon = u.T, v.T, p.T, epsilon.T
    return u, v, p


####
def run_simulation():
    """
    Run the simulation.

    Args:
        params: Dictionary containing the simulation parameters.
    """

    t_start = 0
    t_end = 3
    dt_io = 1
    folder_path = "results/"
    dx = 0.1
    xmin = 0
    xmax = 40
    ymin = 0
    ymax = 10
    cfl = 0.4

    # Initialize the 2D grid
    X, Y = np.arange(xmin, xmax, dx), np.arange(ymin, ymax, dx)
    xx, yy = np.meshgrid(X, Y)
    NX, NY = len(X), len(Y)

    # potential mask
    epsilon = np.zeros((NX, NY))

    lx = NX * dx
    ly = NY * dx
    kx = np.array([(2 * np.pi * i / lx) for i in range(0, (int(NX / 2) - 1))])
    kx = np.append(kx, np.array([(2 * np.pi * (NX - i) / lx)
                   for i in range(int(NX / 2) - 1, NX)]))
    ky = np.array([(np.pi * (i + 1) / ly) for i in range(0, NY)])
    KX, KY = np.meshgrid(kx, ky)
    K = KX ** 2 + KY ** 2

    ###########################################
    u_init = np.ones((NX, NY))
    v_init = np.zeros_like(u_init)
    p_init = np.ones_like(u_init)

    u = copy.deepcopy(u_init)
    v = copy.deepcopy(v_init)
    p = copy.deepcopy(p_init)

    bc = {'x': 'own', 'y': 'free-slip'}
    if bc['y'] == 'no-slip':
        u[0, :] = 0
        u[-1, :] = 0
        v[0, :] = 0
        v[-1, :] = 0

    # Run simulation
    start_time = time.time()
    ts = [t_start]
    t_last = t_start

    # Object in the wind tunnel
    r = ((xx - lx / 4) ** 2 + (yy - ly / 2) ** 2) ** 0.5
    theta = np.arctan2(yy - ly / 2, xx - lx / 4)

    R = 0.5

    # save initial condition
    ip_op.save_step(folder_path, xx, yy, p, u, v, ts[-1])
    # progress_bar = tqdm(total=(t_end - t_start))
    while ts[-1] < t_end:
        uold = copy.deepcopy(u)
        vold = copy.deepcopy(v)
        pold = copy.deepcopy(p)

        # Uold = boundary_conditions(Uold)

        dt = 0.015  # compute_dt(Uold, dx, cfl)
        if ts[-1] + dt >= t_end:
            dt = t_end - ts[-1]
        if dt == 0:
            break
        ts.append(ts[-1] + dt)

        # Update the solution using the indices computed above
        u, v, p = update_step(uold, vold, pold, dt, NX, NY, dx, epsilon, K, bc)

        if ts[-1] - t_last >= dt_io or ts[-1] == t_end:
            ip_op.save_step(folder_path, xx, yy, p, u, v, ts[-1])
            t_last = ts[-1]
            print(t_last)

    end_time = time.time()
    return end_time - start_time


if __name__ == "__main__":
    time = run_simulation()
