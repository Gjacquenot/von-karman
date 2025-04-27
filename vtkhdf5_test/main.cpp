#include <vtkSmartPointer.h>
#include <vtkHDFReader.h>
#include <vtkProbeFilter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>

#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformation.h>

#include <array>
#include <iostream>
#include <vector>
#include <string>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.vtkhdf" << std::endl;
        return EXIT_FAILURE;
    }

    std::string hdfFile = argv[1];

    // Step 1: Define probe points
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0.0, 0.0, 0.0);
    points->InsertNextPoint(1.0, 0.0, 0.0);

    auto probeInput = vtkSmartPointer<vtkPolyData>::New();
    probeInput->SetPoints(points);
    int numPoints = points->GetNumberOfPoints();

    // Step 2: Set up HDF reader
    auto reader = vtkSmartPointer<vtkHDFReader>::New();
    reader->SetFileName(hdfFile.c_str());
    reader->UpdateInformation(); // Load metadata to access time steps
    vtkInformation* info = reader->GetOutputInformation(0);

    int numTimeSteps = info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    std::vector<double> timeSteps(numTimeSteps);
    info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0]);

    // Container: allTimeStepData[time][point][component]
    std::vector<std::vector<std::array<double, 3>>> allTimeStepData(
        numTimeSteps, std::vector<std::array<double, 3>>(numPoints));

    // Step 3: Loop over all time steps
    for (int t = 0; t < numTimeSteps; ++t) {
        double currentTime = timeSteps[t];
        reader->UpdateTimeStep(currentTime);
        reader->Update();

        auto probe = vtkSmartPointer<vtkProbeFilter>::New();
        probe->SetInputData(probeInput);
        probe->SetSourceConnection(reader->GetOutputPort());
        probe->Update();

        vtkDataArray* vecArray = probe->GetOutput()->GetPointData()->GetVectors("velocity");
        if (!vecArray) {
            std::cerr << "No velocity vector data found at timestep " << currentTime << std::endl;
            continue;
        }

        for (int i = 0; i < numPoints; ++i) {
            double vec[3];
            vecArray->GetTuple(i, vec);
            allTimeStepData[t][i] = {vec[0], vec[1], vec[2]};
        }
    }

    // Step 4: Print collected vector values
    for (int t = 0; t < numTimeSteps; ++t) {
        std::cout << "Time " << timeSteps[t] << ":\n";
        for (int i = 0; i < numPoints; ++i) {
            const auto& v = allTimeStepData[t][i];
            std::cout << "  Point " << i << ": (" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
        }
    }

    return EXIT_SUCCESS;
}
