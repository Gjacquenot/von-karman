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
#include <iostream>

int main()
{
    // Assume this is your input multiblock with CellData vectors
    vtkSmartPointer<vtkMultiBlockDataSet> inputMultiBlock = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    // Create a simple unstructured grid and add it to the multiblock dataset
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0.0, 0.0, 0.0);
    points->InsertNextPoint(1.0, 0.0, 0.0);
    points->InsertNextPoint(0.0, 1.0, 0.0);
    grid->SetPoints(points);

    // Define a triangular cell using the points
    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
    triangle->GetPointIds()->SetId(0, 0);
    triangle->GetPointIds()->SetId(1, 1);
    triangle->GetPointIds()->SetId(2, 2);

    grid->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());

    vtkSmartPointer<vtkDoubleArray> vectors = vtkSmartPointer<vtkDoubleArray>::New();
    vectors->SetNumberOfComponents(3);
    vectors->SetName("Velocity");
    vectors->InsertNextTuple3(1.0, 0.0, 0.0);
    grid->GetCellData()->SetVectors(vectors);

    inputMultiBlock->SetBlock(0, grid);

    // First pass: collect valid blocks
    std::vector<vtkDataSet*> validBlocks;
    unsigned int numBlocks = inputMultiBlock->GetNumberOfBlocks();

    for (unsigned int i = 0; i < numBlocks; ++i)
    {
        vtkDataSet* ds = vtkDataSet::SafeDownCast(inputMultiBlock->GetBlock(i));
        if (!ds || !ds->GetCellData()->GetVectors()) continue;
        validBlocks.push_back(ds);
    }

    if (validBlocks.empty())
    {
        std::cerr << "No valid blocks found with cell data vectors!" << std::endl;
        return -1;
    }

    // Second pass: convert CellData to PointData
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

    // Set up the interpolated velocity field
    vtkSmartPointer<vtkCompositeInterpolatedVelocityField> velocityField =
        vtkSmartPointer<vtkCompositeInterpolatedVelocityField>::New();

    // Use point data
    velocityField->Initialize(pointDataMultiBlock, 1);  // 1 = use point data

    // Add each block from pointDataMultiBlock to velocityField
    for (unsigned int i = 0; i < pointDataMultiBlock->GetNumberOfBlocks(); ++i)
    {
        vtkDataSet* pointDataSet = vtkDataSet::SafeDownCast(pointDataMultiBlock->GetBlock(i));
        if (pointDataSet)
        {
            velocityField->AddDataSet(pointDataSet);
        }
        else
        {
            std::cerr << "Block " << i << " is not a valid vtkDataSet!" << std::endl;
        }
    }

    // Set the vector array name (make sure it exists)
    velocityField->SelectVectors(vtkDataObject::FIELD_ASSOCIATION_POINTS, "Velocity");

    // Setup locator strategy
    vtkSmartPointer<vtkCellLocator> treeLoc = vtkSmartPointer<vtkCellLocator>::New();
    treeLoc->SetDataSet(grid);  // Set the dataset for the locator
    treeLoc->BuildLocator();    // Build the locator

    vtkSmartPointer<vtkCellLocatorStrategy> strategy = vtkSmartPointer<vtkCellLocatorStrategy>::New();
    strategy->SetCellLocator(treeLoc);
    velocityField->SetFindCellStrategy(strategy);

    // Sample position inside the bounds (adjust as needed)
    double x[3] = {0.1, 0.1, 0.0};  // pick a point inside your domain
    double u[3] = {0.0, 0.0, 0.0};
    vtkIndent indent;
    velocityField->PrintSelf(std::cout, indent);
    // Evaluate the vector field at the given point
    int result = velocityField->FunctionValues(x, u);
    if (result)
    {
        std::cout << "Interpolated velocity at [" << x[0] << ", " << x[1] << ", " << x[2] << "] is: "
                  << u[0] << ", " << u[1] << ", " << u[2] << std::endl;
    }
    else
    {
        std::cerr << "Failed to evaluate velocity field at point!" << std::endl;
    }

    return 0;
}
