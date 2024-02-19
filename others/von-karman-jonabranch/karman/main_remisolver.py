import sys
import importlib.util
#import numpy as np
#import matplotlib.pyplot as plt
#import h5py
#import copy
#from tqdm import tqdm
from worker_simulation import run_simulation
#from helper_setup import *
#from helper_setup import *


def create_simulation_description(params, t_numba):
    initial_condition_name = params["initial_conditions"]
    boundary_conditions_name = params["boundary_conditions"]
    object_condition_name = params["object_condition"]
    t_start = params["t_start"]
    t_end = params["t_end"]
    dt = params["dt"]
    dt_io = params["dt_io"]
    filename_save = params["filename_save"]
    dx = params["dx"]
    xmin = params["xmin"]
    xmax = params["xmax"]
    ymin = params["ymin"]
    ymax = params["ymax"]
    gamma = params["gamma"]
    cfl = params["cfl"]
    
    
    # Create a text file and write setup information
    with open(filename_save.split(".")[0] + "_description.txt", "w") as f:
        f.write("Simulation Description:\n")
        f.write(f"Initial Conditions: {initial_condition_name}\n")
        f.write(f"Boundary Conditions: {boundary_conditions_name}\n")
        f.write(f"Object: {object_condition_name}\n")
        f.write(f"Start Time: {t_start}\n")
        f.write(f"End Time: {t_end}\n")
        f.write(f"Time Step Size: {dt} (0 iff cfl condition)\n")
        f.write(f"Output Time Interval: {dt_io}\n")
        f.write(f"Spatial Step Size: {dx}\n")
        f.write(f"X-axis Range: {xmin} to {xmax}\n")
        f.write(f"Y-axis Range: {ymin} to {ymax}\n")
        f.write(f"Gamma: {gamma}\n")
        f.write(f"CFL: {cfl}\n")
        f.write(f"Filename Save: {filename_save}\n")

        f.write("Runtime: {t_numba}\n")
    return None

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python main_simulation.py setup_script")
        sys.exit(1)

    setup_script = sys.argv[1]

    # Import simulation parameters from the setup script
    spec = importlib.util.spec_from_file_location("setup_simulation", setup_script)
    setup_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(setup_module)

    # Extract simulation parameters
    params = {
        "initial_conditions": setup_module.initial_conditions,
        "boundary_conditions": setup_module.boundary_conditions,
        "object_condition": setup_module.object_condition,
        "t_start": setup_module.t_start,
        "t_end": setup_module.t_end,
        "dt": setup_module.dt,
        "dt_io": setup_module.dt_io,
        "filename_save": setup_module.filename_save,
        "dx": setup_module.dx,
        "xmin": setup_module.xmin,
        "xmax": setup_module.xmax,
        "ymin": setup_module.ymin,
        "ymax": setup_module.ymax,
        "gamma": setup_module.gamma,
        "cfl": setup_module.cfl
    }

    # Run the simulation
    t_numba = run_simulation(params)
    create_simulation_description(params, t_numba)