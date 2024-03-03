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
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from defaults import *
from plot_functions import *


def animate(folder_path: str, folder_object, obstacle, save: bool = False):
    # Create a figure
    fig, ax = plt.subplots()

    # Read data
    num_frames = get_num_frames(folder_path)

    # set data
    times, data_blocks_u, data_blocks_v, data_blocks_w, data_blocks_p = set_data(
        folder_path)

    FPS = int(num_frames / duration) + 1  # +1 to ensure positivity
    interval = 1000 / FPS  # interval between frames in milliseconds

    data_blocks = set_datablocks(
        w_on,
        data_blocks_u,
        data_blocks_v,
        data_blocks_w)

    # Print some information
    print("Length of data: ", num_frames)
    print("Number of frames: ", num_frames)
    print("Interval: ", interval)
    print("Expected duration: ", duration)

    # Create data
    nx, ny = set_nx_ny(data_blocks)
    X, Y = set_axis(dx, dy, LX, LY, nx, ny)

    Z = data_blocks[0]
    Z_MAX, Z_MIN = set_Z_max_min(data_blocks, w_on)

    color = get_color(w_on)

    plot_args = get_plot_args(dx, LX, LY, Z_MIN, Z_MAX, color)

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
    if obstacle != "":
        get_obstacle(obstacle, nx, ax, folder_object)

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
folder_results = 'output/results/'
folder_object = 'config/'

# Read data
Re, NX, NY, LX, LY, L, U, nu, dx, dy, dt, T, obstacle, w_on, animation_on = readSetupFromHDF5(
    folder_results + '../setup.h5')

if animation_on == False:
    print("Animation is off. Exiting...")
    sys.exit()

# data = read_data_object(folder_object + obstacle + '.txt')
# if obstacle == 'circle':
#     x0, y0, radius = data
#     print(x0, y0, radius)
# elif obstacle == 'rectangle':
#     x0, y0, width, height = data
#     print(x0, y0, width, height)
# elif obstacle == "mountain":
#     x0, y0, lamb, h = data
#     print(x0, y0, lamb, h)
# elif obstacle == "airfoil":
#     a, b, c, d, e, lamb, x0, y0 = data
#     print(a, b, c, d, e, lamb, x0, y0)


animate(folder_results, folder_object, obstacle)
