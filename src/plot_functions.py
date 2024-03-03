import numpy as np
from matplotlib import colormaps as cm
from matplotlib.colors import ListedColormap
import os
from read_data import read_data_object
from matplotlib.patches import Circle, Rectangle
from defaults import *


def get_color(w_on):
    redblue = cm['RdBu_r']
    only_red = ListedColormap(redblue(np.linspace(0.5, 1, 128)))
    color = redblue if w_on else only_red
    return color


def set_nx_ny(data_blocks):
    nx = data_blocks.shape[2]
    ny = data_blocks.shape[1]
    return nx, ny


def set_axis(dx, dy, LX, LY, nx, ny):
    X = np.linspace(dx / 2 * (1 + cols_removed_beginning),
                    LX - dx / 2 * (1 + cols_removed_end), nx)
    Y = np.linspace(dy / 2, LY - dy / 2, ny)
    X, Y = np.meshgrid(X, Y, indexing='xy')
    return X, Y


def get_plot_args(dx, LX, LY, Z_MIN, Z_MAX, color):
    plot_args = {
        'cmap': color,
        'extent': [
            dx * cols_removed_beginning,
            LX - dx * cols_removed_end,
            0,
            LY],
        'vmin': Z_MIN,
        'vmax': Z_MAX,
        'interpolation': 'spline16'}
    return plot_args


def set_datablocks(w_on, data_blocks_u, data_blocks_v, data_blocks_w):
    if w_on:
        data_blocks = data_blocks_w
    else:
        data_blocks = np.sqrt(data_blocks_u**2 + data_blocks_v**2)
    return data_blocks


def set_Z_max_min(data_blocks, w_on):
    Z_MAX = np.max(data_blocks)
    Z_MIN = np.min(data_blocks)
    if w_on:
        Z_MAX = max(Z_MAX, -Z_MIN) / 4
        Z_MIN = -Z_MAX
    return Z_MAX, Z_MIN


def get_obstacle(obstacle, nx, ax, folder_object):
    data = read_data_object(folder_object + obstacle + '.txt')
    if obstacle == 'circle':
        x0, y0, radius = data
        obstacle = Circle((x0, y0), radius, color=color_object, fill=True)
        ax.add_artist(obstacle)
    elif obstacle == 'rectangle':
        x0, y0, width, height = data
        obstacle = Rectangle(
            (x0 - width / 2, y0 - height / 2), width, height, color=color_object, fill=True)
        ax.add_artist(obstacle)
    elif obstacle == 'mountain':
        # mountain is the 1D function f(x) = y0 - sqrt(lambda ^ 2 (x - x0) ^ 2
        # + h ^ 2)
        x0, y0, lamb, h = data
        XX = np.linspace(
            x0 -
            np.sqrt(
                (y0 ** 2 -
                 h ** 2) /
                lamb ** 2),
            x0 +
            np.sqrt(
                (y0 ** 2 -
                 h ** 2) /
                lamb ** 2),
            5 * nx)
        YY1 = y0 + np.sqrt(y0 - np.sqrt(lamb**2 * (XX - x0)**2 + h**2))
        YY2 = 0
        obstacle = ax.fill_between(XX, YY1, YY2, color=color_object)
    elif obstacle == 'airfoil':
        a, b, c, d, e, lamb, x0, y0 = data
        # airfoil is the 1D function f(x) = y0 + a * sqrt(x - x0) + b * (x
        # - x0) + c * (x - x0)^2 + d * (x - x0)^3 + e * (x - x0)^4
        XX = np.linspace(x0, x0 + 1, 5 * nx)
        YY1 = y0 + a * np.sqrt(XX - x0) + b * (XX - x0) + c * \
            (XX - x0)**2 + d * (XX - x0)**3 + e * (XX - x0)**4
        YY2 = y0 - lamb * (a * np.sqrt(XX - x0) + b * (XX - x0) +
                           c * (XX - x0)**2 + d * (XX - x0)**3 + e * (XX - x0)**4)
        obstacle = ax.fill_between(XX, YY1, YY2, color=color_object)
