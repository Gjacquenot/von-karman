// g++ -std=c++11 demo.cpp -o demo -lhdf5_cpp -lhdf5
// g++ -std=c++11 demo_vtkhdf.cpp -o demo -isystem /usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5
// Copyright: Guillaume Jacquenot
// Licence: MIT
#include <iostream>
#include <vector>
#include <cmath>
#include "H5Cpp.h"

const H5std_string FILE_NAME("vtk_output_temporal_cpp.hdf");
const H5std_string DATASET_NAME("VTKHDF");
const H5std_string POINT_DATA_NAME("PointData");
const H5std_string STEPS_NAME("Steps");
const H5std_string SCALAR_DATASET_NAME("TemporalScalars");
const H5std_string VECTOR_DATASET_NAME("TemporalVectors");

int main() {
    // Define grid dimensions
    const int nx = 40, ny = 50, nz = 60;
    const int num_time_steps = 51;

    // Define metadata for the VTK HDF5 file
    const double origin[3] = {0.0, 0.0, 0.0};
    const double spacing[3] = {0.1, 0.1, 0.1};
    const int extent[6] = {0, nx - 1, 0, ny - 1, 0, nz - 1};
    const double direction[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    // Define temporal data
    std::vector<double> time_steps(num_time_steps);
    for (int t = 0; t < num_time_steps; ++t) {
        time_steps[t] = static_cast<double>(t) * 10.0 / (num_time_steps - 1);
    }

    // Create the HDF5 file
    H5::H5File file(FILE_NAME, H5F_ACC_TRUNC);

    // Create the VTKHDF group
    H5::Group vtkhdf_group = file.createGroup(DATASET_NAME);
    hsize_t dims[1] = {2};
    H5::DataSpace attr_dataspace(1, dims);
    H5::Attribute version_attr = vtkhdf_group.createAttribute("Version", H5::PredType::NATIVE_INT, attr_dataspace);
    int version[2] = {2, 0};
    version_attr.write(H5::PredType::NATIVE_INT, version);

    hsize_t dims_type[1] = {1};
    H5::DataSpace attr_dataspace_type(1, dims_type);
    H5::Attribute type_attr = vtkhdf_group.createAttribute("Type", H5::StrType(H5::PredType::C_S1, 9), attr_dataspace_type);
    const char* type = "ImageData";
    type_attr.write(H5::StrType(H5::PredType::C_S1, 9), type);

    dims[0] = 3;
    H5::DataSpace attr_dataspace_3d(1, dims);
    H5::Attribute origin_attr = vtkhdf_group.createAttribute("Origin", H5::PredType::NATIVE_DOUBLE, attr_dataspace_3d);
    origin_attr.write(H5::PredType::NATIVE_DOUBLE, origin);

    H5::Attribute spacing_attr = vtkhdf_group.createAttribute("Spacing", H5::PredType::NATIVE_DOUBLE, attr_dataspace_3d);
    spacing_attr.write(H5::PredType::NATIVE_DOUBLE, spacing);

    dims[0] = 6;
    H5::DataSpace attr_dataspace_6d(1, dims);
    H5::Attribute extent_attr = vtkhdf_group.createAttribute("WholeExtent", H5::PredType::NATIVE_INT, attr_dataspace_6d);
    extent_attr.write(H5::PredType::NATIVE_INT, extent);

    dims[0] = 9;
    H5::DataSpace attr_dataspace_9d(1, dims);
    H5::Attribute direction_attr = vtkhdf_group.createAttribute("Direction", H5::PredType::NATIVE_DOUBLE, attr_dataspace_9d);
    direction_attr.write(H5::PredType::NATIVE_DOUBLE, direction);

    // Create a PointData group
    H5::Group point_data_group = vtkhdf_group.createGroup(POINT_DATA_NAME);
    dims[0] = 1;
    H5::DataSpace point_data_space(1, dims);
    H5::Attribute scalar_attr = point_data_group.createAttribute("Scalars", H5::StrType(H5::PredType::C_S1, 14), point_data_space);
    scalar_attr.write(H5::StrType(H5::PredType::C_S1, 14), SCALAR_DATASET_NAME);

    H5::Attribute vector_attr = point_data_group.createAttribute("Vectors", H5::StrType(H5::PredType::C_S1, 14), point_data_space);
    vector_attr.write(H5::StrType(H5::PredType::C_S1, 14), VECTOR_DATASET_NAME);

    // Create the Steps group for temporal data
    H5::Group steps_group = vtkhdf_group.createGroup(STEPS_NAME);
    dims[0] = 1;
    H5::Attribute nsteps_attr = steps_group.createAttribute("NSteps", H5::PredType::NATIVE_INT, attr_dataspace);
    nsteps_attr.write(H5::PredType::NATIVE_INT, &num_time_steps);

    hsize_t time_dims[1] = {num_time_steps};
    H5::DataSpace time_dataspace(1, time_dims);
    H5::DataSet time_dataset = steps_group.createDataSet("Values", H5::PredType::NATIVE_DOUBLE, time_dataspace);
    time_dataset.write(time_steps.data(), H5::PredType::NATIVE_DOUBLE);

    // Initialize scalar and vector datasets with temporal dimensions
    hsize_t scalar_dims[4] = {static_cast<hsize_t>(num_time_steps), static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx)};
    H5::DataSpace scalar_dataspace(4, scalar_dims);
    H5::DataSet scalar_dataset = point_data_group.createDataSet(SCALAR_DATASET_NAME, H5::PredType::NATIVE_DOUBLE, scalar_dataspace);

    hsize_t vector_dims[5] = {static_cast<hsize_t>(num_time_steps), static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx), 3};
    H5::DataSpace vector_dataspace(5, vector_dims);
    H5::DataSet vector_dataset = point_data_group.createDataSet(VECTOR_DATASET_NAME, H5::PredType::NATIVE_DOUBLE, vector_dataspace);

    // Write temporal data
    std::vector<double> scalar_data(nx * ny * nz);
    std::vector<double> vector_data(nx * ny * nz * 3);
    for (int t_idx = 0; t_idx < num_time_steps; ++t_idx) {
        double t = time_steps[t_idx];
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int idx = k * ny * nx + j * nx + i;
                    scalar_data[idx] = sin(t + idx);
                    vector_data[3 * idx] = cos(t + 3 * idx);
                    vector_data[3 * idx + 1] = cos(t + 3 * idx + 1);
                    vector_data[3 * idx + 2] = cos(t + 3 * idx + 2);
                }
            }
        }

        // Define the hyperslab for the current time step
        hsize_t start[4] = {static_cast<hsize_t>(t_idx), 0, 0, 0};
        hsize_t count[4] = {1, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx)};
        scalar_dataspace.selectHyperslab(H5S_SELECT_SET, count, start);

        hsize_t vector_start[5] = {static_cast<hsize_t>(t_idx), 0, 0, 0, 0};
        hsize_t vector_count[5] = {1, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx), 3};
        vector_dataspace.selectHyperslab(H5S_SELECT_SET, vector_count, vector_start);

        // Define memory space for scalar data
        hsize_t scalar_mem_dims[3] = {static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx)};
        H5::DataSpace scalar_memspace(3, scalar_mem_dims);

        // Define memory space for vector data
        hsize_t vector_mem_dims[4] = {static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx), 3};
        H5::DataSpace vector_memspace(4, vector_mem_dims);

        // Write data to the datasets
        scalar_dataset.write(scalar_data.data(), H5::PredType::NATIVE_DOUBLE, scalar_memspace, scalar_dataspace);
        vector_dataset.write(vector_data.data(), H5::PredType::NATIVE_DOUBLE, vector_memspace, vector_dataspace);
    }

    std::cout << "VTKHDF file '" << FILE_NAME << "' with temporal point data has been created." << std::endl;
    return 0;
}