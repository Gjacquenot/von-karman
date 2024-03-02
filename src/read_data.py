import numpy as np
import h5py


# def read_data_file(file_path):
#     """Read data from file. The file is of the form:
#     header_0
#     x_00 x_01 ... x_0m
#     x_10 x_11 ... x_1m
#     ...
#     x_n0 x_n1 ... x_nm

#     header_1
#     x_00 x_01 ... x_0m
#     x_10 x_11 ... x_1m
#     ...
#     x_n0 x_n1 ... x_nm

#     ...

#     header_k
#     x_00 x_01 ... x_0m
#     x_10 x_11 ... x_1m
#     ...
#     x_n0 x_n1 ... x_nm
#     """
#     headers = []
#     data_blocks = []

#     with open(file_path, 'r') as file:
#         lines = file.readlines()

#     i = 0
#     while i < len(lines):
#         line = lines[i].strip()
#         if line:
#             header = float(line)
#             headers.append(header)

#             i += 1
#             block_lines = []
#             while i < len(lines) and lines[i].strip():
#                 aux = list(map(float, lines[i].split()))
#                 block_lines.append(aux)
#                 i += 1
#             data_block = np.array(block_lines)
#             data_blocks.append(data_block)
#         i += 1

#     headers_array = np.array(headers)
#     data_blocks_array = np.array(data_blocks)

#     return headers_array, data_blocks_array


def readSetupFromHDF5(hdf5_file):
    with h5py.File(hdf5_file, 'r') as file:
        # Accessing the first element and converting to scalar
        Re = file['Re'][()].item()
        NX = int(file['NX'][()].item())
        NY = int(file['NY'][()].item())
        LX = file['LX'][()].item()
        LY = file['LY'][()].item()
        L = file['L'][()].item()
        U = file['U'][()].item()
        nu = file['nu'][()].item()
        dx = file['dx'][()].item()
        dy = file['dy'][()].item()
        dt = file['dt'][()].item()
        T = file['T'][()].item()
        w_on = bool(file['w_on'][()].item())
        animation_on = bool(file['animation_on'][()].item())
    return Re, NX, NY, LX, LY, L, U, nu, dx, dy, dt, T, w_on, animation_on


def readSolutionFromHDF5(hdf5_file):
    with h5py.File(hdf5_file, 'r') as file:
        # NX = int(file['NX'][()].item())
        # NY = int(file['NY'][()].item())
        u = file['u'][:]
        v = file['v'][:]
        w = file['w'][:]
        p = file['p'][:]
        t = file['t'][()].item()
    return t, u, v, w, p
