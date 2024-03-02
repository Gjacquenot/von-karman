import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from matplotlib import colormaps as cm
from matplotlib.ticker import LinearLocator
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.axes import Axes
from matplotlib.text import Text
from typing import cast
import sys
import time
from read_data import *
from matplotlib.patches import Circle, Rectangle


# DEFAULTS
X_LABEL = 'x'
Y_LABEL = 'y'
Z_LABEL = 'z'
FONTSIZE_TIME = 10
y_pos_text = 1.03
y_pos_title = y_pos_text + 0.14
colorbar_args = {'shrink': 0.5, 'aspect': 10, 'location': 'left'}
duration = 10  # duration of the animation in seconds
cols_removed_beginning = 0
cols_removed_end = 20


def animate(folder_path: str, save: bool = False):
    # Create a figure
    fig, ax = plt.subplots()

    # Read data
    num_frames = len([f for f in os.listdir(folder_path)
                      if os.path.isfile(os.path.join(folder_path, f))])

    # Read data
    Re, NX, NY, LX, LY, L, U, nu, dx, dy, dt, T, w_on, animation_on = readSetupFromHDF5(
        folder_path + '../setup.h5')
    times = []
    data_blocks_u = []
    data_blocks_v = []
    data_blocks_w = []
    data_blocks_p = []
    for i in range(num_frames):
        file_path = folder_path + 'sol_' + str(i) + '.h5'
        t, u, v, w, p = readSolutionFromHDF5(file_path)
        times.append(t)
        # remove ghost cells
        data_blocks_u.append(u[1:-1, 1:-1])
        data_blocks_v.append(v[1:-1, 1:-1])
        data_blocks_w.append(w[1:-1, 1:-1])
        data_blocks_p.append(p[1:-1, 1:-1])

    times = np.array(times)
    data_blocks_u = np.array(data_blocks_u)
    data_blocks_v = np.array(data_blocks_v)
    data_blocks_w = np.array(data_blocks_w)
    data_blocks_p = np.array(data_blocks_p)

    FPS = int(num_frames / duration) + 1  # +1 to ensure positivity
    interval = 1000 / FPS  # interval between frames in milliseconds

    # data_blocks = np.sqrt(data_blocks_u**2 + data_blocks_v**2)
    data_blocks = data_blocks_w

    # Print some information
    print("Length of data: ", num_frames)
    print("Number of frames: ", num_frames)
    print("Interval: ", interval)
    print("Expected duration: ", duration)

    # Create data
    nx_old = data_blocks.shape[1]
    ny_old = data_blocks.shape[2]

    # delete the first cols_removed_beginning and last cols_removed_end columns
    data_blocks = np.delete(data_blocks, range(cols_removed_beginning), axis=1)
    data_blocks = np.delete(data_blocks, range(-cols_removed_end, 0), axis=1)

    # Create data
    nx = data_blocks.shape[1]
    ny = data_blocks.shape[2]

    dx = LX / nx_old
    dy = LY / ny_old

    X = np.linspace(dx / 2 * (1 + cols_removed_beginning),
                    LX - dx / 2 * (1 + cols_removed_end), nx)
    Y = np.linspace(dy / 2, LY - dy / 2, ny)
    X, Y = np.meshgrid(X, Y, indexing='xy')

    # transpose the data and reverse
    old_data_blocks = data_blocks
    data_blocks = np.zeros((len(old_data_blocks), ny, nx))
    for i in range(len(old_data_blocks)):
        data_blocks[i] = old_data_blocks[i].T[::-1]

    Z = data_blocks[0]
    Z_MAX = np.max(data_blocks)
    Z_MIN = np.min(data_blocks)
    if w_on:
        # the division by 2 is to make the colors more visible
        Z_MAX = max(Z_MAX, -Z_MIN) / 2
        Z_MIN = -Z_MAX
    # concentrate more strong colors (red and blue) on lower values
    color = cm['RdBu_r'] if w_on else cm['inferno']
    plot_args = {
        'cmap': color,
        'extent': [
            dx * cols_removed_beginning,
            LX - dx * cols_removed_end,
            0,
            LY],
        'vmin': Z_MIN,
        'vmax': Z_MAX,
        # 'interpolation': 'none'}
        'interpolation': 'spline16'}
    # 'interpolation': 'bilinear'}
    text_args = {'x': 0.5, 'y': y_pos_text, 's': '', 'transform': ax.transAxes,
                 'fontsize': FONTSIZE_TIME, 'horizontalalignment': 'center'}
    time_text = ax.text(**text_args)

    # set title (omega or sqrt(u^2 + v^2)
    ax.set_title(
        r'$\omega$' if w_on else r'$\sqrt{u^2 + v^2}$',
        y=y_pos_title)

    # Axes
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)

    # plot a circle
    if object != "":
        if object == 'circle':
            if w_on:
                obstacle = Circle((x0, y0), radius, color='black', fill=True)
            else:
                obstacle = Circle((x0, y0), radius, color='w', fill=False)
            ax.add_artist(obstacle)
        elif object == 'rectangle':
            if w_on:
                obstacle = Rectangle(
                    (x0 - width / 2, y0 - height / 2), width, height, color='black', fill=True)
            else:
                obstacle = Rectangle(
                    (x0 - width / 2, y0 - height / 2), width, height, color='w', fill=False)
            ax.add_artist(obstacle)
        elif object == 'mountain':
            # mountain is the 1D function f(x) = y0 - sqrt(lambda^2 (x - x0)^2
            # + h)
            if w_on:
                obstacle = ax.plot(X[0],
                                   y0 - np.sqrt(lamb**2 * (X[0] - x0)**2 + h),
                                   color='black')
            else:
                obstacle = ax.plot(X[0],
                                   y0 - np.sqrt(lamb**2 * (X[0] - x0)**2 + h),
                                   color='w')
        elif object == 'airfoil':
            # airfoil is the 1D function f(x) = y0 + a * sqrt(x - x0) + b * (x
            # - x0) + c * (x - x0)^2 + d * (x - x0)^3 + e * (x - x0)^4
            XX = np.linspace(x0, x0 + 1, 5 * nx)
            if w_on:
                obstacle = ax.plot(XX, y0 +
                                   a * np.sqrt(XX - x0) +
                                   b * (XX - x0) +
                                   c * (XX - x0)**2 +
                                   d * (XX - x0)**3 +
                                   e * (XX - x0)**4,
                                   color='black') + ax.plot(XX, y0 -
                                                            lamb * (a * np.sqrt(XX - x0) +
                                                                    b * (XX - x0) +
                                                                    c * (XX - x0)**2 +
                                                                    d * (XX - x0)**3 +
                                                                    e * (XX - x0)**4),
                                                            color='black')

            else:
                obstacle = ax.plot(XX, y0 +
                                   a * np.sqrt(XX - x0) +
                                   b * (XX - x0) +
                                   c * (XX - x0)**2 +
                                   d * (XX - x0)**3 +
                                   e * (XX - x0)**4,
                                   color='w') + ax.plot(XX, y0 -
                                                        lamb * (a * np.sqrt(XX - x0) +
                                                                b * (XX - x0) +
                                                                c * (XX - x0)**2 +
                                                                d * (XX - x0)**3 +
                                                                e * (XX - x0)**4),
                                                        color='w')
    # Create plot
    plot = [ax.imshow(Z, **plot_args)]

    # Add a color bar which maps values to colors.
    fig.colorbar(plot[0], **colorbar_args)

    # Animation update function
    def init():
        return plot[0], time_text

    def update(frame):
        # Clear the previous frame
        # ax.clear()
        plot[0].remove()

        real_frame = int(frame)

        # Update the arrays
        Z = data_blocks[real_frame]

        # Update the time text
        time_text.set_text('t = %.3f' % times[real_frame])

        # Update the plot
        plot[0] = ax.imshow(Z, **plot_args)

        return plot[0], time_text

    # Create the animation
    ani = FuncAnimation(fig, update, frames=num_frames,
                        interval=interval, blit=False, init_func=init)

    end_time = time.time()
    print("Total time for animating: ", int(
        UNIT_TIME * (end_time - start_time)), LABEL_TIME)

    if save:
        # Save the animation
        ani.save('animation.mp4', writer='ffmpeg', fps=FPS)
    else:
        # Show the animation
        plt.show()


# count time
UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()

# print python executable instruction
# Join the arguments using a space as a separator
arguments = ' '.join(sys.argv[1:])
print("python src/animation.py " + arguments)

try:
    object = sys.argv[1]
except IndexError:
    object = ""

if object == 'circle':
    try:
        x0 = float(sys.argv[2])
        y0 = float(sys.argv[3])
        radius = float(sys.argv[4])
    except IndexError:
        x0 = np.nan
        y0 = np.nan
        radius = np.nan
elif object == 'rectangle':
    try:
        x0 = float(sys.argv[2])
        y0 = float(sys.argv[3])
        width = float(sys.argv[4])
        height = float(sys.argv[5])
    except IndexError:
        x0 = np.nan
        y0 = np.nan
        width = np.nan
        height = np.nan
elif object == "mountain":
    try:
        x0 = float(sys.argv[2])
        y0 = float(sys.argv[3])
        h = float(sys.argv[4])
        lamb = float(sys.argv[5])
    except IndexError:
        x0 = np.nan
        y0 = np.nan
        h = np.nan
        lamb = np.nan
elif object == "airfoil":
    try:
        a = float(sys.argv[2])
        b = float(sys.argv[3])
        c = float(sys.argv[4])
        d = float(sys.argv[5])
        e = float(sys.argv[6])
        lamb = float(sys.argv[7])
        x0 = float(sys.argv[8])
        y0 = float(sys.argv[9])
    except IndexError:
        a = np.nan
        b = np.nan
        c = np.nan
        d = np.nan
        e = np.nan
        lamb = np.nan
        x0 = np.nan
        y0 = np.nan

folder_name = 'output/results/'

animate(folder_name)
