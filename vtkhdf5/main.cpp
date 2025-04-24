#include "vtkhdf5.hpp"
#include <vector>
#include <cmath>

int main() {
    const int nx = 40, ny = 50, nz = 2; // 2D grid (nz = 1)

    // Initialize writer
    VTKHDFWriter writer("vtk_output_temporal_cpp.hdf", nx, ny, nz);

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
        writer.writeTimestep(time, scalar_data, vector_data);
    }

    return 0;
}