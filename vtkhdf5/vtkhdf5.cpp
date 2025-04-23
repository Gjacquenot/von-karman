#include "vtkhdf5.hpp"
#include <iostream>

VTKHDFWriter::VTKHDFWriter(const std::string &filename, int nx, int ny, int nz)
    : file(filename, H5F_ACC_TRUNC), nx(nx), ny(ny), nz(nz), current_timestep(0) {
    const double origin[3] = {0.0, 0.0, 0.0};
    const double spacing[3] = {0.1, 0.1, 0.1};
    const int extent[6] = {0, nx - 1, 0, ny - 1, 0, nz - 1};
    const double direction[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    // Create the VTKHDF group
    H5::Group vtkhdf_group = file.createGroup("VTKHDF");

    // Add the Type attribute
    hsize_t dims_type[1] = {1};
    H5::DataSpace attr_dataspace_type(1, dims_type);
    H5::Attribute type_attr = vtkhdf_group.createAttribute("Type", H5::StrType(H5::PredType::C_S1, 9), attr_dataspace_type);
    const char *type = "ImageData";
    type_attr.write(H5::StrType(H5::PredType::C_S1, 9), type);

    // Metadata attributes
    hsize_t dims[1] = {2};
    H5::DataSpace attr_dataspace(1, dims);
    int version[2] = {2, 0};
    vtkhdf_group.createAttribute("Version", H5::PredType::NATIVE_INT, attr_dataspace)
        .write(H5::PredType::NATIVE_INT, version);

    dims[0] = 3;
    H5::DataSpace attr_dataspace_3d(1, dims);
    vtkhdf_group.createAttribute("Origin", H5::PredType::NATIVE_DOUBLE, attr_dataspace_3d)
        .write(H5::PredType::NATIVE_DOUBLE, origin);
    vtkhdf_group.createAttribute("Spacing", H5::PredType::NATIVE_DOUBLE, attr_dataspace_3d)
        .write(H5::PredType::NATIVE_DOUBLE, spacing);

    dims[0] = 6;
    H5::DataSpace attr_dataspace_6d(1, dims);
    vtkhdf_group.createAttribute("WholeExtent", H5::PredType::NATIVE_INT, attr_dataspace_6d)
        .write(H5::PredType::NATIVE_INT, extent);

    dims[0] = 9;
    H5::DataSpace attr_dataspace_9d(1, dims);
    vtkhdf_group.createAttribute("Direction", H5::PredType::NATIVE_DOUBLE, attr_dataspace_9d)
        .write(H5::PredType::NATIVE_DOUBLE, direction);

    // Create PointData group
    H5::Group point_data_group = vtkhdf_group.createGroup("PointData");
    dims[0] = 1;
    H5::DataSpace point_data_space(1, dims);
    point_data_group.createAttribute("Scalars", H5::StrType(H5::PredType::C_S1, 14), point_data_space)
        .write(H5::StrType(H5::PredType::C_S1, 14), "TemporalScalars");
    point_data_group.createAttribute("Vectors", H5::StrType(H5::PredType::C_S1, 14), point_data_space)
        .write(H5::StrType(H5::PredType::C_S1, 14), "TemporalVectors");

    // Create Steps group
    H5::Group steps_group = vtkhdf_group.createGroup("Steps");
    // Create time dataset with unlimited dimension
    hsize_t time_dims[1] = {1};
    hsize_t max_time_dims[1] = {H5S_UNLIMITED};
    H5::DataSpace time_dataspace(1, time_dims, max_time_dims);
    H5::DSetCreatPropList time_prop;
    time_prop.setChunk(1, time_dims);
    time_dataset = steps_group.createDataSet("Values", H5::PredType::NATIVE_DOUBLE, time_dataspace, time_prop);

    // Create scalar and vector datasets with unlimited time dimension
    hsize_t scalar_dims[4] = {1, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx)};
    hsize_t max_scalar_dims[4] = {H5S_UNLIMITED, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx)};
    scalar_dataspace = H5::DataSpace(4, scalar_dims, max_scalar_dims);
    H5::DSetCreatPropList scalar_prop;
    scalar_prop.setChunk(4, scalar_dims);
    scalar_dataset = point_data_group.createDataSet("TemporalScalars", H5::PredType::NATIVE_DOUBLE, scalar_dataspace, scalar_prop);

    hsize_t vector_dims[5] = {1, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx), 3};
    hsize_t max_vector_dims[5] = {H5S_UNLIMITED, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx), 3};
    vector_dataspace = H5::DataSpace(5, vector_dims, max_vector_dims);
    H5::DSetCreatPropList vector_prop;
    vector_prop.setChunk(5, vector_dims);
    vector_dataset = point_data_group.createDataSet("TemporalVectors", H5::PredType::NATIVE_DOUBLE, vector_dataspace, vector_prop);
}

VTKHDFWriter::~VTKHDFWriter() {
    try {
        // Write the NSteps attribute to the Steps group
        H5::Group steps_group = file.openGroup("VTKHDF/Steps");
        hsize_t dims[1] = {1};
        H5::DataSpace attr_dataspace(1, dims);
        H5::Attribute nsteps_attr = steps_group.createAttribute("NSteps", H5::PredType::NATIVE_INT, attr_dataspace);
        int nsteps = static_cast<int>(current_timestep);
        nsteps_attr.write(H5::PredType::NATIVE_INT, &nsteps);
    } catch (H5::Exception &e) {
        std::cerr << "Error writing NSteps attribute: " << e.getDetailMsg() << std::endl;
    }
}

void VTKHDFWriter::writeTimestep(double time, const std::vector<double> &scalar_data, const std::vector<double> &vector_data) {
    // Extend datasets to accommodate the new timestep
    hsize_t new_time_dims[1] = {current_timestep + 1};
    time_dataset.extend(new_time_dims);

    hsize_t new_scalar_dims[4] = {current_timestep + 1, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx)};
    scalar_dataset.extend(new_scalar_dims);

    hsize_t new_vector_dims[5] = {current_timestep + 1, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx), 3};
    vector_dataset.extend(new_vector_dims);

    // Write time value
    hsize_t time_start[1] = {current_timestep};
    H5::DataSpace time_memspace(1, new_time_dims);
    time_dataset.write(&time, H5::PredType::NATIVE_DOUBLE, time_memspace, time_dataset.getSpace());

    // Write scalar data
    hsize_t scalar_start[4] = {current_timestep, 0, 0, 0}; // Start at the current timestep
    hsize_t scalar_count[4] = {1, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx)}; // 3D grid dimensions
    H5::DataSpace scalar_memspace(4, scalar_count); // Memory space for the 3D grid
    scalar_dataspace = scalar_dataset.getSpace(); // Update dataspace after extending the dataset
    scalar_dataspace.selectHyperslab(H5S_SELECT_SET, scalar_count, scalar_start); // Select hyperslab for the current timestep
    scalar_dataset.write(scalar_data.data(), H5::PredType::NATIVE_DOUBLE, scalar_memspace, scalar_dataspace);

    // Write vector data
    hsize_t vector_start[5] = {current_timestep, 0, 0, 0, 0}; // Start at the current timestep
    hsize_t vector_count[5] = {1, static_cast<hsize_t>(nz), static_cast<hsize_t>(ny), static_cast<hsize_t>(nx), 3}; // 3D grid dimensions + vector components
    H5::DataSpace vector_memspace(5, vector_count); // Memory space for the 3D grid with vector components
    vector_dataspace = vector_dataset.getSpace(); // Update dataspace after extending the dataset
    vector_dataspace.selectHyperslab(H5S_SELECT_SET, vector_count, vector_start); // Select hyperslab for the current timestep
    vector_dataset.write(vector_data.data(), H5::PredType::NATIVE_DOUBLE, vector_memspace, vector_dataspace);

    // Increment timestep counter
    current_timestep++;
}