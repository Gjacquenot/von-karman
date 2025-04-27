cmake_minimum_required(VERSION 3.16)
project(ProbeHDFVelocity)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find VTK (requires version >= 9.2 for vtkHDFReader)
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    IOHDF
    IOLegacy
    FiltersCore
    FiltersGeneral
    FiltersSources
    FiltersHybrid
    FiltersGeometry
    FiltersPoints
    ImagingCore
    IOImage
    InteractionStyle
    RenderingCore
    RenderingFreeType
    RenderingOpenGL2
)

# Include VTK headers
include(${VTK_USE_FILE})

# Add executable
add_executable(ProbeHDFVelocity main.cpp)

# Link to VTK libraries
target_link_libraries(ProbeHDFVelocity PRIVATE ${VTK_LIBRARIES})
