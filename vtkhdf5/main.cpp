// g++ -std=c++11 main.cpp vtkhdf5.cpp -o demo -isystem /usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5
// Copyright: Guillaume Jacquenot
// Licence: MIT
#include "vtkhdf5.hpp"
#include <array>
#include <cmath>
#include <vector>

int demo_2d() {
    const int nx = 40, ny = 50;

    const std::array<int, 2> dimensions({nx, ny});
    VTKHDFWriter2D writer("vtk_output_temporal_cpp_2d.hdf", dimensions);
    // Simulate and write timesteps
    for (int t_idx = 0; t_idx < 51; ++t_idx) {
        double time = static_cast<double>(t_idx) * 10.0 / 50;

        // Generate scalar and vector data
        std::vector<double> scalar_data(nx * ny);
        std::vector<double> vector_data(nx * ny * 3);
        const double period = 10.0;
        const double omega = 2.0 * M_PI / period;
        const double kx = omega * omega / 9.81;
        const double ky = 0;
        int idx = 0;
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int ij = j * nx + i;
                scalar_data[idx] = sin(-omega * time + kx * i + ky * j);
                vector_data[3 * idx] = cos(time + 3 * ij);
                vector_data[3 * idx + 1] = cos(time + 3 * ij + 1);
                vector_data[3 * idx + 2] = cos(time + 3 * ij + 2);
                ++idx;
        }
        }

        // Write timestep
        writer.write_timestep_vectors(time, scalar_data, vector_data);
    }

    return 0;
}

int demo_3d() {
    const int nx = 40, ny = 50, nz = 2; // 2D grid (nz = 1)
    const std::array<int, 3> nx_ny_nz({nx, ny, nz});
    // Initialize writer
    VTKHDFWriter3D writer("vtk_output_temporal_cpp_3d.hdf", nx_ny_nz);
    // Simulate and write timesteps
    for (int t_idx = 0; t_idx < 51; ++t_idx) {
        double time = static_cast<double>(t_idx) * 10.0 / 50;

        // Generate scalar and vector data
        std::vector<double> scalar_data(nx * ny * nz);
        std::vector<double> vector_data(nx * ny * nz * 3);
        int idx = 0;
        for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int ij = 0*(k * ny * nx) + j * nx + i;
                scalar_data[idx] = sin(time + ij);
                vector_data[3 * idx] = cos(time + 3 * ij);
                vector_data[3 * idx + 1] = cos(time + 3 * ij + 1);
                vector_data[3 * idx + 2] = cos(time + 3 * ij + 2);
                ++idx;
            }
        }
        }

        // Write timestep
        writer.write_timestep_vectors(time, scalar_data, vector_data);
    }

    return 0;
}


int main() {
    demo_2d();
    demo_3d();

    return 0;
}