#include <vtkSmartPointer.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDataSet.h>
#include <vtkCellDataToPointData.h>
#include <vtkCompositeInterpolatedVelocityField.h>
#include <vtkCellLocatorStrategy.h>
#include <vtkCellLocator.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkDataObject.h>
#include <vtkUnstructuredGrid.h>
#include <vtkNew.h>
#include <vtkTriangle.h>
#include <vtkProbeFilter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkIndent.h>
#include <iostream>

int main()
{
    vtkSmartPointer<vtkMultiBlockDataSet> inputMultiBlock = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    // Create a sample unstructured grid
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0.0, 0.0, 0.0);
    points->InsertNextPoint(1.0, 0.0, 0.0);
    points->InsertNextPoint(0.0, 1.0, 0.0);
    grid->SetPoints(points);

    // Triangle cell
    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
    triangle->GetPointIds()->SetId(0, 0);
    triangle->GetPointIds()->SetId(1, 1);
    triangle->GetPointIds()->SetId(2, 2);
    grid->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());

    // Velocity (vector) data
    vtkSmartPointer<vtkDoubleArray> vectors = vtkSmartPointer<vtkDoubleArray>::New();
    vectors->SetNumberOfComponents(3);
    vectors->SetName("Velocity");
    vectors->InsertNextTuple3(1.0, 0.0, 0.0);
    vectors->InsertNextTuple3(1.0, 0.0, 0.0);
    vectors->InsertNextTuple3(1.0, 0.0, 0.0);
    grid->GetCellData()->SetVectors(vectors);

    // Pressure (scalar) data
    vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
    scalars->SetName("Pressure");
    scalars->InsertNextValue(101.0);
    scalars->InsertNextValue(102.0);
    scalars->InsertNextValue(103.0);
    grid->GetCellData()->SetScalars(scalars);

    inputMultiBlock->SetBlock(0, grid);

    // Filter valid blocks
    std::vector<vtkDataSet*> validBlocks;
    unsigned int numBlocks = inputMultiBlock->GetNumberOfBlocks();
    for (unsigned int i = 0; i < numBlocks; ++i)
    {
        vtkDataSet* ds = vtkDataSet::SafeDownCast(inputMultiBlock->GetBlock(i));
        if (ds && ds->GetCellData()->GetVectors())
            validBlocks.push_back(ds);
    }

    if (validBlocks.empty())
    {
        std::cerr << "No valid blocks with vector data!" << std::endl;
        return -1;
    }

    // Convert to point data
    vtkSmartPointer<vtkMultiBlockDataSet> pointDataMultiBlock = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    pointDataMultiBlock->SetNumberOfBlocks(validBlocks.size());

    for (size_t i = 0; i < validBlocks.size(); ++i)
    {
        vtkSmartPointer<vtkCellDataToPointData> converter = vtkSmartPointer<vtkCellDataToPointData>::New();
        converter->SetInputData(validBlocks[i]);
        converter->PassCellDataOff();
        converter->Update();

        vtkDataSet* pointDataSet = converter->GetOutput();
        pointDataMultiBlock->SetBlock(i, pointDataSet);
    }

    // ----------- Velocity Field (Vector) -----------
    vtkSmartPointer<vtkCompositeInterpolatedVelocityField> velocityField =
        vtkSmartPointer<vtkCompositeInterpolatedVelocityField>::New();
    velocityField->Initialize(pointDataMultiBlock, 1); // 1 = use point data
    velocityField->SelectVectors(vtkDataObject::FIELD_ASSOCIATION_POINTS, "Velocity");

    for (unsigned int i = 0; i < pointDataMultiBlock->GetNumberOfBlocks(); ++i)
    {
        vtkDataSet* ds = vtkDataSet::SafeDownCast(pointDataMultiBlock->GetBlock(i));
        if (ds) velocityField->AddDataSet(ds);
    }

    // Locator strategy for velocity field
    vtkSmartPointer<vtkCellLocator> velLocator = vtkSmartPointer<vtkCellLocator>::New();
    velLocator->SetDataSet(grid);
    velLocator->BuildLocator();
    vtkSmartPointer<vtkCellLocatorStrategy> velStrategy = vtkSmartPointer<vtkCellLocatorStrategy>::New();
    velStrategy->SetCellLocator(velLocator);
    velocityField->SetFindCellStrategy(velStrategy);

    // ----------- Pressure Field (Scalar) using vtkProbeFilter -----------
    vtkSmartPointer<vtkCellDataToPointData> scalarConverter = vtkSmartPointer<vtkCellDataToPointData>::New();
    scalarConverter->SetInputData(grid);
    scalarConverter->PassCellDataOff();
    scalarConverter->Update();

    vtkSmartPointer<vtkPolyData> probePoint = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> probePts = vtkSmartPointer<vtkPoints>::New();
    probePts->InsertNextPoint(0.1, 0.1, 0.0);
    probePoint->SetPoints(probePts);

    vtkSmartPointer<vtkProbeFilter> probe = vtkSmartPointer<vtkProbeFilter>::New();
    probe->SetSourceData(scalarConverter->GetOutput());
    probe->SetInputData(probePoint);
    probe->Update();

    vtkDataArray* pressureValues = probe->GetOutput()->GetPointData()->GetArray("Pressure");

    // ----------- Test evaluation -----------
    double x[3] = {0.1, 0.1, 0.0};
    double u[3] = {0.0, 0.0, 0.0};

    int valid = velocityField->FunctionValues(x, u);
    if (valid)
    {
        std::cout << "Velocity @ [0.1, 0.1, 0.0] = (" << u[0] << ", " << u[1] << ", " << u[2] << ")\n";
    }
    else
    {
        std::cerr << "Failed to evaluate velocity field.\n";
    }

    if (pressureValues)
    {
        std::cout << "Pressure @ [0.1, 0.1, 0.0] = " << pressureValues->GetComponent(0, 0) << "\n";
    }
    else
    {
        std::cerr << "Failed to probe pressure.\n";
    }

    return 0;
}
