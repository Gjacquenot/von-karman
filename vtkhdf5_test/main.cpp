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


int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.vtkhdf" << std::endl;
        return EXIT_FAILURE;
    }

    std::string hdfFile = argv[1];

    // Step 1: Read vtkHDF image data
    auto reader = vtkSmartPointer<vtkHDFReader>::New();
    reader->SetFileName(hdfFile.c_str());
    reader->UpdateInformation();

    // Get time steps
    vtkInformation* info = reader->GetOutputInformation(0);
    int numTimeSteps = info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    std::cout << "Number of time steps = " << numTimeSteps << std::endl;
    std::vector<double> timeSteps(numTimeSteps);
    if (numTimeSteps > 0)
    info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0]);
    for (int t = 0; t < numTimeSteps; ++t) {
        std::cout << "t("<<t <<") "<< timeSteps[t] << std::endl;
    }

    // Read first timestep to get image metadata
    if (numTimeSteps > 0)
        reader->UpdateTimeStep(timeSteps[0]);
    else
        reader->Update();

    vtkImageData* image = vtkImageData::SafeDownCast(reader->GetOutput());

    double bounds[6];
    image->GetBounds(bounds);

    std::cout << "📦 Bounding box (world coordinates):\n";
    std::cout << "  X: [" << bounds[0] << ", " << bounds[1] << "]\n";
    std::cout << "  Y: [" << bounds[2] << ", " << bounds[3] << "]\n";
    std::cout << "  Z: [" << bounds[4] << ", " << bounds[5] << "]\n";

    // Optional: define some probe points
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0.0, 0.0, 0.0);
    points->InsertNextPoint(1.0, 0.0, 0.0);

    auto probeInput = vtkSmartPointer<vtkPolyData>::New();
    probeInput->SetPoints(points);
    int numPoints = points->GetNumberOfPoints();

    std::vector<std::vector<std::array<double, 3>>> allTimeStepData(
        numTimeSteps, std::vector<std::array<double, 3>>(numPoints));

    for (int t = 0; t < numTimeSteps; ++t) {
        reader->UpdateTimeStep(timeSteps[t]);
        reader->Update();

        auto probe = vtkSmartPointer<vtkProbeFilter>::New();
        probe->SetInputData(probeInput);
        probe->SetSourceConnection(reader->GetOutputPort());
        probe->Update();

        vtkDataArray* vecArray = probe->GetOutput()->GetPointData()->GetVectors("TemporalVectors");
        if (!vecArray) {
            std::cerr << "No 'velocity' vector data at t=" << timeSteps[t] << std::endl;
            continue;
        }

        for (int i = 0; i < numPoints; ++i) {
            double v[3];
            vecArray->GetTuple(i, v);
            allTimeStepData[t][i] = {v[0], v[1], v[2]};
        }
    }

    // Display interpolated data
    for (int t = 0; t < numTimeSteps; ++t) {
        std::cout << "🕒 Time " << timeSteps[t] << ":\n";
        for (int i = 0; i < numPoints; ++i) {
            const auto& v = allTimeStepData[t][i];
            std::cout << "  Point " << i << ": (" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
        }
    }

    // Set to a specific time step
    size_t desiredIndex =10;
    double time = timeSteps[desiredIndex];
    reader->UpdateTimeStep(time);
    reader->Update();

    // Get output and cast safely
    vtkImageData* image2 = vtkImageData::SafeDownCast(reader->GetOutput());
    if (!image2) {
        std::cerr << "Error: Output is not vtkImageData!" << std::endl;
        return EXIT_FAILURE;
    }

    // Get attribute from point data
    // vtkDataArray* velocityArray = image2->GetPointData()->GetVectors("velocity");
    vtkDataArray* velocityArray = image2->GetPointData()->GetVectors("TemporalVectors");
    if (!velocityArray) {
        std::cerr << "No 'velocity' array found!" << std::endl;
        return EXIT_FAILURE;
    }

    vtkIdType numTuples = velocityArray->GetNumberOfTuples();
    int numComponents = velocityArray->GetNumberOfComponents();

    std::cout << "Dumping velocity vectors (" << numTuples << " points):\n";

    Timer timer;
    double v[3];
    for (vtkIdType i = 0; i < numTuples; ++i) {
        velocityArray->GetTuple(i, v);
        // std::cout << "Point " << i << ": (" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
    }
    timer.stopAndDisplay("Total computation time");

    return EXIT_SUCCESS;
}
