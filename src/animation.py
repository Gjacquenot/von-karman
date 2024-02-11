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
from read_data import read_data_file


# DEFAULTS
X_LABEL = 'x'
Y_LABEL = 'y'
Z_LABEL = 'z'
FONTSIZE_TIME = 10
y_pos_text = 1.03
y_pos_title = y_pos_text + 0.14
color = cm['inferno']
colorbar_args = {'shrink': 0.5, 'aspect': 10, 'location': 'left'}
FPS = 50
interval = 1000 / FPS  # interval between frames in milliseconds
duration = 40  # duration of the animation in seconds


def animate(save=False):
    # Create a figure
    fig, ax = plt.subplots()

    # Read data
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Replace with the actual file path
    file_path = script_dir + '/../data/solution.txt'
    headers, data_blocks = read_data_file(file_path)

    # in practice it takes more time to animate, so we divide the duration
    # by 2 to keep the seconds more or less the same

    # frames to be animated
    num_frames = int(duration * 1000.0 / interval)
    if num_frames > len(data_blocks):
        num_frames = len(data_blocks)

    # Speed of the animation
    speed = len(data_blocks) / num_frames

    print("Number of frames: ", num_frames)
    print("Speed: ", speed)
    print("Interval: ", interval)
    print("Expected duration: ", duration)
    # print("len(data_blocks): ", len(data_blocks))
    # print("len(headers): ", len(headers))
    # print("len(data_blocks_freq): ", len(data_blocks_freq))

    # Create data
    nx = data_blocks.shape[1]
    ny = data_blocks.shape[2]

    dx = Lx / nx
    dy = Ly / ny

    X = np.linspace(dx / 2, Lx - dx / 2, nx)
    Y = np.linspace(dy / 2, Ly - dy / 2, ny)
    X, Y = np.meshgrid(X, Y, indexing='ij')
    Z = data_blocks[0]
    Z_MAX = np.max(data_blocks)
    Z_MIN = 0
    plot_args = {'cmap': color}
    text_args = {'x': 0.5, 'y': y_pos_text, 's': '', 'transform': ax.transAxes,
                 'fontsize': FONTSIZE_TIME, 'horizontalalignment': 'center'}
    time_text = ax.text(**text_args)

    # Axes
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)

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

        real_frame = int(frame * speed)

        # Update the arrays
        Z = data_blocks[real_frame]

        # Update the time text
        time_text.set_text('t = %.2f' % headers[real_frame])

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
try:
    Lx = float(sys.argv[1])
    Ly = float(sys.argv[2])
except IndexError:
    Lx = np.nan
    Ly = np.nan

animate()
