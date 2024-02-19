import numpy as np
import matplotlib.pyplot as plt
import h5py
import copy
import tqdm

#Related plotting
def contourplot(U, X, Y, title=""):
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    titles = [r'$\rho$', r'$\rho u$', r'$\rho v$', r'$\rho E$']
    for i in range(4):
        ax = axs[i // 2, i % 2]
        im = ax.imshow(U[:, :, i].T, extent=[X.min(), X.max(), Y.min(), Y.max()], aspect='auto', origin='lower', cmap='coolwarm')
        ax.set_title(titles[i])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im, ax=ax, label=titles[i])

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.suptitle(title)
    plt.tight_layout()
    plt.show()




def linesplot(U, X, Y, title=""):
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    titles = [r'Constant x at X=0', r'Constant y at Y=0', r'X = Y', r'X = -Y']

    # Plot constant x slice at X=0
    ax = axs[0, 0]
    x_idx = 0
    slice_idx = np.argmin(np.abs(X - x_idx))
    ax.plot(Y, U[slice_idx, :, :], label = [r'$\rho$', r'$\rho u$', r'$\rho v$', r'$\rho E$'])
    ax.set_title(titles[0])
    ax.set_xlabel('Y')
    ax.set_ylabel('Value')
    ax.legend()
    ax.grid(True)

    # Plot constant y slice at Y=0
    ax = axs[0, 1]
    y_idx = 0
    slice_idx = np.argmin(np.abs(Y - y_idx))
    ax.plot(X, U[:, slice_idx, :], label = [r'$\rho$', r'$\rho u$', r'$\rho v$', r'$\rho E$'])
    ax.set_title(titles[1])
    ax.set_xlabel('X')
    ax.set_ylabel('Value')
    ax.legend()
    ax.grid(True)
    if len(X) == len(Y):
        # Plot x=y
        ax = axs[1, 0]
        
        ids = range(len(X))
        
        ax.plot(X, U[ids, ids, :], label = [r'$\rho$', r'$\rho u$', r'$\rho v$', r'$\rho E$'])
        ax.set_title(titles[2])
        ax.set_xlabel('X = Y')
        ax.set_ylabel('Value')
        ax.legend()
        ax.grid(True)

        # Plot x=-y
        ax = axs[1, 1]
    
        slice_idx = np.argmin(np.abs(Y - y_idx))
        ax.plot(X, U[ids, range(U.shape[0]-1, -1, -1) , :], label = [r'$\rho$', r'$\rho u$', r'$\rho v$', r'$\rho E$'])
        ax.set_title(titles[3])
        ax.set_xlabel('X = Y')
        ax.set_ylabel('Value')
        ax.legend()
        ax.grid(True)

    plt.suptitle(title)
    plt.tight_layout()
    plt.show()


def plot_flow_arrows(U, X, Y):
    # Extract velocity components from U
    rho_u = U[:, :, 1]
    rho_v = U[:, :, 2]

    # Create meshgrid for coordinates
    X_mesh, Y_mesh = np.meshgrid(X, Y)

    # Plot arrows
    plt.figure(figsize=(8, 6))
    plt.quiver(X_mesh, Y_mesh, rho_u, rho_v, scale=20)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Flow Arrows')
    plt.show()



def get_conservation(U0, U):
    """
    Calculate conservation metrics based on initial and final states.

    Args:
        U0 (ndarray): Initial state of the system.
        U (ndarray): Final state of the system.

    Returns:
        tuple: Tuple containing conservation metrics for mass, momentum, and energy.
    """
    rhos = np.sum(U[:, :, 0] - U0[:, :, 0]) / np.sum(U0[:, :, 0])
    rhou = np.sum(U[:, :, 1] - U0[:, :, 1]) / np.sum(U0[:, :, 1])
    rhov = np.sum(U[:, :, 2] - U0[:, :, 2]) / np.sum(U0[:, :, 2])
    rhoE = np.sum(U[:, :, 3] - U0[:, :, 3]) / np.sum(U0[:, :, 3])
    return rhos, rhou, rhov, rhoE

def plot_conservation(data):
    """
    Plot conservation metrics for mass, momentum, and energy.

    Args:
        data (dict): Dictionary containing simulation data with time steps as keys.
    """
    conservation = []
    for key in sorted(data.keys()):
        conservation.append(get_conservation(data[0.], data[key]))

    conservation = np.array(conservation)
    fig = plt.figure(figsize=(10, 8))
    plt.plot(sorted(data.keys()), conservation[:, 0], label=r'mean $\rho $')
    plt.plot(sorted(data.keys()), conservation[:, 1], label=r'mean $\rho u $')
    plt.plot(sorted(data.keys()), conservation[:, 2], label=r'mean $\rho v$')
    plt.plot(sorted(data.keys()), conservation[:, 3], label=r'mean $\rho E $')
    plt.xlabel("Time")
    plt.ylabel("Conservation (final_mass - initial_mass) / initial_mass")
    plt.legend()
    plt.grid(True)
    plt.title("Conservation of mass, momentum, and energy")
    plt.show()
