import numpy as np
import os
import shutil
from matplotlib import pyplot
from src import projection_method, helper_setup, ip_op, integrate


def PyDNS(file_path):
    t_start, t_end, dt, dt_io, dx, dy, xmin, xmax, ymin, ymax, cfl, folder_path = helper_setup.read_simulation_setup(
        file_path)
    # free storage
    if os.path.exists(folder_path):
        shutil.rmtree(folder_path)
    os.makedirs(folder_path)
    dt = 0.005
    # variable declarations
    x = np.arange(xmin, xmax, dx)
    nx = len(x)
    y = np.arange(ymin, ymax, dy)
    ny = len(y)
    lx = xmax - xmin
    ly = ymax - ymin
    xx, yy = np.meshgrid(x, y)

    # That needs to be fixed

    F = 0
    dt = 0.015

    # physical variables
    rho = 1
    nu = 1 / 150

    # boundary conditions
    bc = {'x': 'periodic', 'y': 'free-slip'}

    # initial conditions
    u = np.zeros((ny, nx))
    utemp = np.zeros((ny, 3))

    v = np.zeros((ny, nx))
    vtemp = np.zeros((ny, 3))

    if bc['y'] == 'no-slip':
        u[0, :] = 0
        u[-1, :] = 0
        v[0, :] = 0
        v[-1, :] = 0

    p = np.zeros((ny, nx))
    ptemp = np.zeros((ny, 3))
    dpdx = np.zeros((ny, nx))
    dpdy = np.zeros((ny, nx))
    epsilon = np.zeros((ny, nx))
    uRHS_conv_diff_p = np.zeros((ny, nx))
    uRHS_conv_diff_pp = np.zeros((ny, nx))
    vRHS_conv_diff_p = np.zeros((ny, nx))
    vRHS_conv_diff_pp = np.zeros((ny, nx))

    # Object in the wind tunnel
    r = ((xx - lx / 4) ** 2 + (yy - ly / 2) ** 2) ** 0.5
    theta = np.arctan2(yy - ly / 2, xx - lx / 4)

    R = -0.5

    for i in range(nx):
        for j in range(ny):
            if r[j, i] <= R:
                epsilon[j, i] = 1

    # pressure_poisson
    nx_sp = nx
    ny_sp = ny

    kx = np.array([(2 * np.pi * i / lx)
                  for i in range(0, (int(nx_sp / 2) - 1))])
    kx = np.append(kx, np.array([(2 * np.pi * (nx_sp - i) / lx)
                   for i in range(int(nx_sp / 2) - 1, nx_sp)]))
    ky = np.array([(np.pi * (i + 1) / ly) for i in range(0, ny_sp)])
    KX, KY = np.meshgrid(kx, ky)
    K = KX ** 2 + KY ** 2

    ip_op.save_step(folder_path, xx, yy, p, u, v, 0)
    t_last = 0

    # rk3
    ts = [0]
    for stepcount in range(1, 3):
        ts.append(ts[-1] + dt)
        u, v, p = projection_method.BC(u, v, p)
        u, v, dpdx, dpdy, uRHS_conv_diff, vRHS_conv_diff = integrate.rk3(u, v, nx, ny, nu, dx, dy, dt, dpdx, dpdy,
                                                                         epsilon, theta, r, R,
                                                                         rho, x, y, xx, yy,
                                                                         nx_sp, ny_sp, K, bc)

        # print("Step=%06i time=%4.6f" % (stepcount, stepcount * dt))

        # if (np.mod(stepcount, saveth_iter) == 0) and (stepcount > save_start):
        #    print("snapshot= %i" % (stepcount / saveth_iter))
        #    ip_op.write_szl_2D(xx, yy, p, u, v, stepcount * dt, int(stepcount / saveth_iter))

        ip_mdot = dy * ((u[0, 0] + u[-1, 0]) / 2 + sum(u[1:-1, 0]))
        op_mdot = dy * ((u[0, -1] + u[-1, -1]) / 2 + sum(u[1:-1, -1]))
        print("mass flow rate ip op diff: %f %f %e" %
              (ip_mdot, op_mdot, op_mdot - ip_mdot))

        uRHS_conv_diff_pp = uRHS_conv_diff_p.copy()
        vRHS_conv_diff_pp = vRHS_conv_diff_p.copy()

        uRHS_conv_diff_p = uRHS_conv_diff.copy()
        vRHS_conv_diff_p = vRHS_conv_diff.copy()

    while ts[-1] < t_end:
        if ts[-1] + 0.015 > t_end:
            ts.append(t_end)
        else:
            ts.append(ts[-1] + dt)
        u, v, p = projection_method.BC(u, v, p)
        ustar, vstar, uRHS_conv_diff, vRHS_conv_diff = projection_method.step1(u, v, nx, ny, nu, x, y, xx, yy, dx, dy,
                                                                               dt, epsilon, R, theta, r,
                                                                               uRHS_conv_diff_p, uRHS_conv_diff_pp,
                                                                               vRHS_conv_diff_p, vRHS_conv_diff_pp,
                                                                               dpdx, dpdy, bc)

        # Step2
        ustarstar, vstarstar = projection_method.step2(
            ustar, vstar, dpdx, dpdy, dt)
        ustarstar, vstarstar, p = projection_method.BC(ustarstar, vstarstar, p)
        # Step3
        p = projection_method.step3(
            ustarstar,
            vstarstar,
            rho,
            epsilon,
            dx,
            dy,
            nx_sp,
            ny_sp,
            K,
            dt,
            bc)

        # Step4
        u, v, dpdx, dpdy = projection_method.step4(
            ustarstar, vstarstar, p, dx, dy, dt, bc)

        if bc['y'] == 'free-slip':
            u[0, :] = u[1, :].copy()
            u[-1, :] = u[-2, :].copy()
            v[0, :] = 0
            v[-1, :] = 0

        # print("Step=%06i time=%4.6f" % (stepcount, stepcount * dt))
        if ts[-1] - t_last >= dt_io or ts[-1] == t_end:
            ip_op.save_step(folder_path, xx, yy, p, u, v, ts[-1])
            t_last = ts[-1]
            print(t_last)

        ip_mdot = dy * ((u[0, 0] + u[-1, 0]) / 2 + sum(u[1:-1, 0]))
        op_mdot = dy * ((u[0, -1] + u[-1, -1]) / 2 + sum(u[1:-1, -1]))
        # print("mass flow rate ip op diff: %f %f %e" % (ip_mdot, op_mdot, op_mdot-ip_mdot))

        uRHS_conv_diff_pp = uRHS_conv_diff_p.copy()
        vRHS_conv_diff_pp = vRHS_conv_diff_p.copy()

        uRHS_conv_diff_p = uRHS_conv_diff.copy()
        vRHS_conv_diff_p = vRHS_conv_diff.copy()

        '''
        fig = pyplot.figure(figsize=(11, 7), dpi=100)
        pyplot.quiver(xx[::3, ::3], yy[::3, ::3], u[::3, ::3], v[::3, ::3])

        fig = pyplot.figure(figsize=(11, 7), dpi=100)
        pyplot.quiver(xx, yy, u, v)

        pyplot.show()
        '''


if __name__ == "__main__":
    import os
    import sys
    from pathlib import Path
    from urllib import request

    if len(sys.argv) != 2:
        print("Usage: python PyDNS.py setup_script.txt")
        sys.exit(1)

    setup_script = sys.argv[1]
    """
    # To download tecio library module
    lib_folder = Path("tecio")
    if os.name == 'nt':
        dll_path = lib_folder / 'libtecio.dll'
        if not dll_path.is_file():
            url = 'https://raw.githubusercontent.com/blacksong/pytecio/master/2017r3_tecio.dll'
            print('Downloading dll from github:', url)
            request.urlretrieve(url, dll_path)

    else:
        dll_path = lib_folder / 'libtecio.so'
        if not dll_path.is_file():
            url = 'https://raw.githubusercontent.com/blacksong/pytecio/master/2017r2_tecio.so'
            print('Downloading dll from github:', url)
            request.urlretrieve(url, dll_path)
    """

    PyDNS(setup_script)
