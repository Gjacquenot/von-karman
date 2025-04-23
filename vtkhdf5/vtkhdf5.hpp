#ifndef VTKHDF_WRITER_HPP
#define VTKHDF_WRITER_HPP

#include <H5Cpp.h>
#include <vector>
#include <string>

class VTKHDFWriter {
public:
    VTKHDFWriter(const std::string &filename, int nx, int ny, int nz);
    ~VTKHDFWriter();

    void writeTimestep(double time, const std::vector<double> &scalar_data, const std::vector<double> &vector_data);

private:
    H5::H5File file;
    H5::DataSet scalar_dataset;
    H5::DataSet vector_dataset;
    H5::DataSet time_dataset;
    H5::DataSpace scalar_dataspace;
    H5::DataSpace vector_dataspace;

    int nx, ny, nz;
    hsize_t current_timestep;
};

#endif // VTKHDF_WRITER_HPP