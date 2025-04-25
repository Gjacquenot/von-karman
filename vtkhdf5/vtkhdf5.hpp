#ifndef VTKHDF_WRITER_HPP
#define VTKHDF_WRITER_HPP

#include <H5Cpp.h>
#include <vector>
#include <string>

class AbstractVTKHDFWriter {
    public:
        AbstractVTKHDFWriter(const std::string& filename);
        virtual ~AbstractVTKHDFWriter();
        void write_timestep_data(double time);
        virtual void write_timestep(double time, const std::vector<double> &scalar_data, const std::vector<double> &vector_data) = 0;
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
        VTKHDFWriter2D(const std::string &filename, int nx, int ny);
        ~VTKHDFWriter2D();
        void write_timestep(double time, const std::vector<double> &scalar_data, const std::vector<double> &vector_data);

    private:
        int nx, ny;
        VTKHDFWriter2D();
};

class VTKHDFWriter3D: public AbstractVTKHDFWriter {
    public:
        VTKHDFWriter3D(const std::string &filename, int nx, int ny, int nz);
        ~VTKHDFWriter3D();
        void write_timestep(double time, const std::vector<double> &scalar_data, const std::vector<double> &vector_data);

    private:
        int nx, ny, nz;
        VTKHDFWriter3D();
};

#endif // VTKHDF_WRITER_HPP