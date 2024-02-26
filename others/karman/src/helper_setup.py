def read_simulation_setup(file_path):
    params = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                key_value_pair = line.split('#')[0].strip().split('=')
                if len(key_value_pair) == 2:
                    key, value = key_value_pair
                    params[key.strip()] = value.strip()
    # Convert string values to appropriate types if necessary
    initial_condition_name = params["initial_conditions"]
    boundary_conditions_name = params["boundary_conditions"]

    t_start= float(params["t_start"])
    t_end = float(params["t_end"])
    dt = float(params["dt"])
    dt_io = float(params["dt_io"])
    dx = float(params["dx"])
    dy = float(params["dy"])
    xmin = float(params["xmin"])
    xmax = float(params["xmax"])
    ymin = float(params["ymin"])
    ymax = float(params["ymax"])
    gamma = float(params["gamma"])
    cfl = float(params["cfl"])
    folder_path = params["folder_path"]
    return t_start, t_end, dt, dt_io, dx, dy, xmin, xmax, ymin, ymax, cfl, folder_path 


