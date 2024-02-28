import numpy as np
from src import derive, pressure_poisson, ibm
import copy


def step1(u, v, nx, ny, nu, x, y, xx, yy, dx, dy, dt, epsilon, R, theta, r, uRHS_conv_diff_p, uRHS_conv_diff_pp,
          vRHS_conv_diff_p, vRHS_conv_diff_pp, dpdx, dpdy, bc):

    # Step1
    # do the x-momentum RHS
    # u rhs: - d(uu)/dx - d(vu)/dy + ν d2(u)
    uRHS_conv_diff = derive.conv_diff_u(u, v, nu, dx, dy, bc)
    # v rhs: - d(uv)/dx - d(vv)/dy + ν d2(v)
    vRHS_conv_diff = derive.conv_diff_v(u, v, nu, dx, dy, bc)

    uRHS = (23 * uRHS_conv_diff - 16 * uRHS_conv_diff_p +
            5 * uRHS_conv_diff_pp) / 12 - dpdx

    vRHS = (23 * vRHS_conv_diff - 16 * vRHS_conv_diff_p +
            5 * vRHS_conv_diff_pp) / 12 - dpdy

    ibm_forcing_u, ibm_forcing_v = ibm.circle(
        uRHS, vRHS, u, v, dt, x, y, xx, yy, r, R, theta, nx, ny, epsilon)
    # uprim = Semilag2(u, v, u, dx, dy, dt, nx, ny)
    # vprim = Semilag2(u, v, v, dx, dy, dt, nx, ny)
    ustar = u + ibm_forcing_u * dt + dt * uRHS
    vstar = v + ibm_forcing_v * dt + dt * vRHS

    if bc['y'] == 'no-slip':
        ustar[0, :] = 0
        ustar[-1, :] = 0
        vstar[0, :] = 0
        vstar[-1, :] = 0
    elif bc['y'] == 'free-slip':
        ustar[0, :] = ustar[1, :].copy()
        ustar[-1, :] = ustar[-2, :].copy()
        v[0, :] = 0
        v[-1, :] = 0

    return ustar, vstar, uRHS_conv_diff, vRHS_conv_diff


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


def Semilag_t(u_in, v_in, q_in, dx, dy, dt, NX, NY):
    """
    1st order semi-Lagrangian advection
    """
    u, v, q = copy.deepcopy(u_in), copy.deepcopy(v_in), copy.deepcopy(q_in)
    u, v, q = u.T, v.T, q.T
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
    return q.T
####


def Semilag2(u, v, q, dx, dy, dt, NX, NY):
    """
    Second order semi-Lagrangian advection
    """
    ADVq = Semilag_t(u, v, q, dx, dy, dt, NX, NY)
    aux = Semilag_t(-u, -v, ADVq, dx, dy, dt, NX, NY)
    ADVq = q + (q - aux) / 2.
    ADVq = Semilag_t(u, v, ADVq, dx, dy, dt, NX, NY)

    return ADVq


def BC(u, v, p):

    u[0, :] = u[1, :].copy()
    u[-1, :] = u[-2, :].copy()
    u[:, 0] = 1
    u[:, -1] = u[:, -2].copy()

    v[:, 0] = 0
    v[:, -1] = v[:, -2].copy()
    v[0, :] = 0
    v[-1, :] = 0

    p[:, 0] = p[:, 1].copy()
    p[:, -1] = 0  # p[:,-2]
    p[0, :] = p[1, :].copy()
    p[-1, :] = p[-2, :].copy()

    return u, v, p


def step2(ustar, vstar, dpdx, dpdy, dt):
    ustarstar = ustar + dpdx * dt
    vstarstar = vstar + dpdy * dt
    return ustarstar, vstarstar


def step3(ustarstar, vstarstar, rho, epsilon, dx, dy, nx_sp, ny_sp, K, dt, bc):
    # next compute the pressure RHS: prhs = div(un)/dt + div( [urhs, vrhs])
    prhs = rho * derive.div((1 - epsilon) * ustarstar,
                            (1 - epsilon) * vstarstar, dx, dy, bc) / dt

    # p, err = pressure_poisson.solve_new(p, dx, dy, prhs)

    p = pressure_poisson.solve_spectral_v2(nx_sp, ny_sp, K, prhs)

    return p


def step4(ustarstar, vstarstar, p, dx, dy, dt, bc):
    # Step4
    # finally compute the true velocities
    # u_{n+1} = uh - dt*dpdx
    dpdx = derive.ddx(p, dx, bc)

    dpdy = derive.ddy(p, dy, bc)

    u = ustarstar - dt * dpdx
    v = vstarstar - dt * dpdy

    return u, v, dpdx, dpdy
