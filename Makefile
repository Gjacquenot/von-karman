# Description: Makefile for the project

# Compiler and flags
CXX := g++
CC := gcc
CXXFLAGS := -Wall -Wextra -Wpedantic -Wshadow -Wformat=2 -Wcast-align -Wconversion -Wsign-conversion -Wnull-dereference -g3 -Ofast -std=c++11
CFLAGS := -Wall -Wextra -pedantic -Ofast -std=c11

# Libraries and paths
HDF5_DIR := /opt/homebrew/Cellar/hdf5/1.14.3
EIGEN_DIR := /opt/homebrew/Cellar/eigen/3.4.0_1
LIBS := -lm -L$(HDF5_DIR)/lib -lhdf5_cpp -lhdf5
INCLUDES := -I$(EIGEN_DIR)/include/ -I$(HDF5_DIR)/include

# Folders
SRC := src
INCLUDE := include
BIN := bin

# Executable name
TARGET := main

# Sources and objects
SOURCES := $(wildcard $(SRC)/*.c) $(wildcard $(SRC)/*.cpp)
OBJECTS := $(patsubst $(SRC)/%.cpp, $(BIN)/%.o, $(filter %.cpp, $(SOURCES))) $(patsubst $(SRC)/%.c, $(BIN)/%.o, $(filter %.c, $(SOURCES)))

# Headers
HEADERS := $(wildcard $(INCLUDE)/*.h) $(wildcard $(INCLUDE)/*.hpp)

.PHONY: all clean

all: $(BIN)/$(TARGET)

$(BIN)/$(TARGET): $(OBJECTS)
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(BIN)/%.o: $(SRC)/%.cpp $(HEADERS)
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BIN)/%.o: $(SRC)/%.c $(HEADERS)
	@$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -f $(BIN)/*.o $(BIN)/$(TARGET)
	@rm -rf $(BIN)
	@rm -rf build

cmake: ## Build with CMake
	@cmake -S . -B build
	@cmake --build build
	@cmake --install build
	@mkdir -p output/results
	# @cmake --build build --target clean