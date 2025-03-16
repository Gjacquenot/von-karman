cmake_minimum_required(VERSION 3.10)
project(von_karman LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable testing
enable_testing()

# Find Eigen3
find_package(Eigen3 REQUIRED)

# Find HDF5 (C++ component)
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Include directories
include_directories(${EIGEN3_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS} include)

# Source files
set(SOURCES src/main.cpp src/NS.cpp)

# Define the executable
add_executable(von_karman ${SOURCES})

# Link libraries
target_link_libraries(von_karman PRIVATE ${HDF5_CXX_LIBRARIES})

# Ensure HDF5 is properly linked
target_include_directories(von_karman PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_directories(von_karman PRIVATE ${HDF5_LIBRARY_DIRS})

# Add a test that runs `von_karman` from the project root
add_test(
  NAME RunVonKarman
  COMMAND ${CMAKE_BINARY_DIR}/von_karman
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
