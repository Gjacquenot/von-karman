import numpy as np
from scipy import interpolate


def circle(uRHS, vRHS, u, v, dt, x, y, xx, yy, r, R, theta, nx, ny, epsilon):
    u_desired = np.zeros_like(u)
    v_desired = np.zeros_like(v)
    interpolate_u = interpolate.interp2d(x, y, u + uRHS * dt, kind='cubic')
    interpolate_v = interpolate.interp2d(x, y, v + vRHS * dt, kind='cubic')

    for i in range(nx):
        for j in range(ny):
            if epsilon[j, i] != 0:
                u_desired[j, i] = -np.sin(2 * np.pi * (r[j, i] ** 2) / (2 * (R ** 2))) * interpolate_u(
                    xx[j, i] + (R - r[j, i]) * np.cos(theta[j, i]),
                    yy[j, i] + (R - r[j, i]) * np.sin(theta[j, i]))
                v_desired[j, i] = -np.sin(2 * np.pi * (r[j, i] ** 2) / (2 * (R ** 2))) * interpolate_v(
                    xx[j, i] + (R - r[j, i]) * np.cos(theta[j, i]),
                    yy[j, i] + (R - r[j, i]) * np.sin(theta[j, i]))

    ibm_forcing_u = epsilon * (-uRHS + (u_desired - u) / dt)
    ibm_forcing_v = epsilon * (-vRHS + (v_desired - v) / dt)

    return ibm_forcing_u, ibm_forcing_v


def ellipse(uRHS, vRHS, u, v, dt, x, y, xx, yy, a, b, nx, ny, epsilon):
    u_desired = np.zeros_like(u)
    v_desired = np.zeros_like(v)
    interpolate_u = interpolate.interp2d(x, y, u + uRHS * dt, kind='cubic')
    interpolate_v = interpolate.interp2d(x, y, v + vRHS * dt, kind='cubic')

    for i in range(nx):
        for j in range(ny):
            # Check if the point (xx, yy) is inside the ellipse
            if (xx[j, i] ** 2 / a ** 2) + (yy[j, i] ** 2 / b ** 2) <= 1:
                # Manipulate the velocity field around the ellipse
                # This is a simple adaptation; more sophisticated models might
                # be required for accurate simulations
                distance = np.sqrt(
                    (xx[j, i] ** 2 / a ** 2) + (yy[j, i] ** 2 / b ** 2))
                u_desired[j, i] = -np.sin(2 * np.pi * (distance ** 2) / 2) * \
                    interpolate_u(xx[j, i], yy[j, i])
                v_desired[j, i] = -np.sin(2 * np.pi * (distance ** 2) / 2) * \
                    interpolate_v(xx[j, i], yy[j, i])

    ibm_forcing_u = epsilon * (-uRHS + (u_desired - u) / dt)
    ibm_forcing_v = epsilon * (-vRHS + (v_desired - v) / dt)

    return ibm_forcing_u, ibm_forcing_v
