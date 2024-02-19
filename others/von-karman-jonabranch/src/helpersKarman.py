"""This module contains helper functions for the von-karman package."""

import numpy as np

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

def set_outer_shell_average(grid, mask):
    """Set the value for the outer shell of the cylinder to the average of the values of the cells outside of the cylinder.

    Args:
        grid (np.ndarray): 2D grid.
        mask (np.ndarray): Mask for the cylinder.
    """
    shell = np.logical_xor(mask, np.roll(mask, 1, axis=0)) | np.logical_xor(mask, np.roll(mask, 1, axis=1))
    shell_values = grid[shell]
    mean_outer_shell = np.mean(shell_values, axis=1)  # Calculate mean along axis 0
    grid[mask] = mean_outer_shell
