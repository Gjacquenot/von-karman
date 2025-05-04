#include <vtkSmartPointer.h>
#include <vtkHDFReader.h>
#include <vtkImageData.h>
#include <vtkProbeFilter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformation.h>

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <stdexcept>
#include <chrono>
#include <iomanip>

class Timer {
public:
    Timer() : start_time(std::chrono::steady_clock::now()) {}

    void reset() {
        start_time = std::chrono::steady_clock::now();
    }

    void stopAndDisplay(const std::string& label = "Elapsed time") const {
        auto end_time = std::chrono::steady_clock::now();
        auto duration = end_time - start_time;

        auto us = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
        double value = us;
        std::string unit = "μs";

        if (value > 1000.0) {
            value /= 1000.0;
            unit = "ms";
        }
        if (value > 1000.0) {
            value /= 1000.0;
            unit = "s";
        }
        if (value > 60.0 && unit == "s") {
            value /= 60.0;
            unit = "min";
        }
        if (value > 60.0 && unit == "min") {
            value /= 60.0;
            unit = "h";
        }

        std::cout << label << ": " << std::fixed << std::setprecision(3) << value << " " << unit << '\n';
    }

private:
    std::chrono::time_point<std::chrono::steady_clock> start_time;
};

class VTKHDFHandler {
public:
    VTKHDFHandler(const std::string& filename) {
        // Initialize reader
        reader = vtkSmartPointer<vtkHDFReader>::New();
        reader->SetFileName(filename.c_str());
        reader->UpdateInformation();

        // Get time steps
        vtkInformation* info = reader->GetOutputInformation(0);
        numTimeSteps = info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        if (numTimeSteps == 0) {
            throw std::runtime_error("No time steps found in the file.");
        }

        timeSteps.resize(numTimeSteps);
        info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0]);

        // Read the first time step to get metadata
        reader->UpdateTimeStep(timeSteps[0]);
        reader->Update();
        vtkImageData* image = vtkImageData::SafeDownCast(reader->GetOutput());
        if (!image) {
            throw std::runtime_error("Failed to read vtkImageData from the file.");
        }

        double bounds[6];
        image->GetBounds(bounds);

        std::cout << "📦 Bounding box (world coordinates):\n";
        std::cout << "  X: [" << bounds[0] << ", " << bounds[1] << "]\n";
        std::cout << "  Y: [" << bounds[2] << ", " << bounds[3] << "]\n";
        std::cout << "  Z: [" << bounds[4] << ", " << bounds[5] << "]\n";
        std::cout << "Number of time steps = " << numTimeSteps << std::endl;
        for (int t = 0; t < numTimeSteps; ++t) {
            std::cout << "t(" << t << ") " << timeSteps[t] << std::endl;
        }

        // List all scalar and vector attributes
        vtkPointData* pointData = image->GetPointData();
        if (!pointData) {
            throw std::runtime_error("No point data found in the file.");
        }

        for (int i = 0; i < pointData->GetNumberOfArrays(); ++i) {
            vtkDataArray* array = pointData->GetArray(i);
            if (!array) continue;

            std::string name = array->GetName();
            if (array->GetNumberOfComponents() == 1) {
                scalarAttributes.push_back(name);
            } else if (array->GetNumberOfComponents() == 3) {
                vectorAttributes.push_back(name);
            }
        }

        std::cout << "Scalar attributes: ";
        for (const auto& attr : scalarAttributes) {
            std::cout << attr << " ";
        }
        std::cout << "\n";

        std::cout << "Vector attributes: ";
        for (const auto& attr : vectorAttributes) {
            std::cout << attr << " ";
        }
        std::cout << "\n";
    }

    std::pair<std::map<std::string, double>, std::map<std::string, std::array<double, 3>>>
    performTemporalSpatialInterpolation(double queryTime, const std::array<double, 3>& queryPoint) {
        // Find the two closest time steps
        int lowerIndex = -1, upperIndex = -1;
        for (int i = 0; i < numTimeSteps - 1; ++i) {
            if (timeSteps[i] <= queryTime && timeSteps[i + 1] >= queryTime) {
                lowerIndex = i;
                upperIndex = i + 1;
                break;
            }
        }

        if (lowerIndex == -1 || upperIndex == -1) {
            throw std::runtime_error("Query time is out of bounds.");
        }

        double t0 = timeSteps[lowerIndex];
        double t1 = timeSteps[upperIndex];
        double alpha = (queryTime - t0) / (t1 - t0);

        // Prepare probe points
        auto queryPoints = vtkSmartPointer<vtkPoints>::New();
        queryPoints->InsertNextPoint(queryPoint.data());

        auto queryPolyData = vtkSmartPointer<vtkPolyData>::New();
        queryPolyData->SetPoints(queryPoints);

        // Probe for t0
        reader->UpdateTimeStep(t0);
        reader->Update();
        auto probeT0 = vtkSmartPointer<vtkProbeFilter>::New();
        probeT0->SetInputData(queryPolyData);
        probeT0->SetSourceConnection(reader->GetOutputPort());
        probeT0->Update();
        vtkPointData* pointDataT0 = probeT0->GetOutput()->GetPointData();

        // Probe for t1
        reader->UpdateTimeStep(t1);
        reader->Update();
        auto probeT1 = vtkSmartPointer<vtkProbeFilter>::New();
        probeT1->SetInputData(queryPolyData);
        probeT1->SetSourceConnection(reader->GetOutputPort());
        probeT1->Update();
        vtkPointData* pointDataT1 = probeT1->GetOutput()->GetPointData();

        // Interpolate scalar and vector attributes
        std::map<std::string, double> interpolatedScalars;
        std::map<std::string, std::array<double, 3>> interpolatedVectors;

        for (const auto& scalarName : scalarAttributes) {
            vtkDataArray* arrayT0 = pointDataT0->GetArray(scalarName.c_str());
            vtkDataArray* arrayT1 = pointDataT1->GetArray(scalarName.c_str());
            if (arrayT0 && arrayT1) {
                double valueT0 = arrayT0->GetTuple1(0);
                double valueT1 = arrayT1->GetTuple1(0);
                double interpolatedValue = (1 - alpha) * valueT0 + alpha * valueT1;
                interpolatedScalars[scalarName] = interpolatedValue;
            }
        }

        for (const auto& vectorName : vectorAttributes) {
            vtkDataArray* arrayT0 = pointDataT0->GetArray(vectorName.c_str());
            vtkDataArray* arrayT1 = pointDataT1->GetArray(vectorName.c_str());
            if (arrayT0 && arrayT1) {
                double v0[3], v1[3];
                arrayT0->GetTuple(0, v0);
                arrayT1->GetTuple(0, v1);
                std::array<double, 3> interpolatedValue = {
                    (1 - alpha) * v0[0] + alpha * v1[0],
                    (1 - alpha) * v0[1] + alpha * v1[1],
                    (1 - alpha) * v0[2] + alpha * v1[2]
                };
                interpolatedVectors[vectorName] = interpolatedValue;
            }
        }

        return {interpolatedScalars, interpolatedVectors};
    }

private:
    vtkSmartPointer<vtkHDFReader> reader;
    std::vector<double> timeSteps;
    std::vector<std::string> scalarAttributes;
    std::vector<std::string> vectorAttributes;
    int numTimeSteps;
};


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.vtkhdf" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        VTKHDFHandler handler(argv[1]);

        Timer timer;
        // Perform temporal-spatial interpolation
        auto result = handler.performTemporalSpatialInterpolation(5.0, {0.5, 0.5, 0.0});

        // Output results
        std::cout << "Interpolated Scalars:\n";
        for (const auto& [name, value] : result.first) {
            std::cout << "  " << name << ": " << value << "\n";
        }

        std::cout << "Interpolated Vectors:\n";
        for (const auto& [name, value] : result.second) {
            std::cout << "  " << name << ": (" << value[0] << ", " << value[1] << ", " << value[2] << ")\n";
        }
        timer.stopAndDisplay("Total computation time");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
