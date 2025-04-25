#ifndef VTKHDF_WRITER_HPP
#define VTKHDF_WRITER_HPP

#include <H5Cpp.h>
#include <array>
#include <vector>
#include <string>

class AbstractVTKHDFWriter {
    public:
        AbstractVTKHDFWriter(const std::string& filename);
        virtual ~AbstractVTKHDFWriter();
        void write_timestep_data(double time);
        void write_timestep_vectors(double time, const std::vector<double> &scalar_data, const std::vector<double> &vector_data);
        virtual void write_timestep(double time, const double * scalar_data, const double * vector_data) = 0;
    protected:
        H5::H5File file;
        H5::Group vtkhdf_group;
        H5::Group point_data_group;
        H5::DataSet scalar_dataset;
        H5::DataSet vector_dataset;
        H5::DataSet time_dataset;
        H5::DataSpace scalar_dataspace;
        H5::DataSpace vector_dataspace;
        hsize_t current_timestep;
};

class VTKHDFWriter2D: public AbstractVTKHDFWriter {
    public:
        VTKHDFWriter2D(const std::string &filename,
            const std::array<int, 2> & nx_ny,
            const std::array<double, 2> &dx_dy = {1.0, 1.0},
            const std::array<double, 2> & ox_oy = {0.0, 0.0});
        ~VTKHDFWriter2D() = default;
        void write_timestep(double time, const double * scalar_data, const double * vector_data);
    private:
        int nx, ny;
        VTKHDFWriter2D();
};

class VTKHDFWriter3D: public AbstractVTKHDFWriter {
    public:
        VTKHDFWriter3D(const std::string &filename,
            const std::array<int, 3> & nx_ny_nz,
            const std::array<double, 3> &dx_dy_dz = {1.0, 1.0, 1.0},
            const std::array<double, 3> &ox_oy_oz = {0.0, 0.0, 0.0});
        ~VTKHDFWriter3D() = default;
        void write_timestep(double time, const double * scalar_data, const double * vector_data);

    private:
        int nx, ny, nz;
        VTKHDFWriter3D();
};

#endif // VTKHDF_WRITER_HPP