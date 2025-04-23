#include "vtkhdf5.hpp"
#include <vector>
#include <cmath>

int main() {
    const int nx = 40, ny = 50, nz = 1; // 2D grid (nz = 1)

    // Initialize writer
    VTKHDFWriter writer("vtk_output_temporal_cpp.hdf", nx, ny, nz);

    // Simulate and write timesteps
    for (int t_idx = 0; t_idx < 51; ++t_idx) {
        double time = static_cast<double>(t_idx) * 10.0 / 50;

        // Generate scalar and vector data
        std::vector<double> scalar_data(nx * ny * nz);
        std::vector<double> vector_data(nx * ny * nz * 3);
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx = j * nx + i;
                scalar_data[idx] = sin(time + idx);
                vector_data[3 * idx] = cos(time + 3 * idx);
                vector_data[3 * idx + 1] = cos(time + 3 * idx + 1);
                vector_data[3 * idx + 2] = cos(time + 3 * idx + 2);
            }
        }

        // Write timestep
        writer.writeTimestep(time, scalar_data, vector_data);
    }

    return 0;
}