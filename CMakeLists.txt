cmake_minimum_required(VERSION 3.18)
project(EikonSolver2D LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to enable CUDA support
option(ENABLE_CUDA "Enable CUDA GPU acceleration" OFF)

if(ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    set(CMAKE_CUDA_ARCHITECTURES 75 80 86)  # RTX 20xx, 30xx, 40xx series
    add_definitions(-DCUDA_ENABLED)
    message(STATUS "CUDA support enabled")
else()
    message(STATUS "CUDA support disabled")
endif()

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Find packages
find_package(Eigen3 REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)

# Add compiler flags for optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
if(ENABLE_CUDA)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --use_fast_math --expt-relaxed-constexpr")
endif()

# Include directories
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(include)

# Create the CPU library
if(ENABLE_CUDA)
    # Include CUDA implementation
    add_library(eikon_solver_2d 
        src/eikon_solver_2d.cpp
        src/eikon_solver_cuda_impl.cu
    )
    # Set CUDA properties
    set_property(TARGET eikon_solver_2d PROPERTY CUDA_STANDARD 17)
    set_property(TARGET eikon_solver_2d PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
    target_link_libraries(eikon_solver_2d Eigen3::Eigen CUDA::cudart)
else()
    # CPU-only version
    add_library(eikon_solver_2d 
        src/eikon_solver_2d.cpp
    )
    target_link_libraries(eikon_solver_2d Eigen3::Eigen)
endif()

# Create HDF5 utility library
add_library(hdf5_utils
    src/hdf5_utils.cpp
)

target_link_libraries(hdf5_utils Eigen3::Eigen ${HDF5_LIBRARIES})

# Create the CPU example executable
add_executable(example_eikon_solver 
    tests/example_eikon_solver.cpp
)

target_link_libraries(example_eikon_solver eikon_solver_2d hdf5_utils Eigen3::Eigen ${HDF5_LIBRARIES})

if(ENABLE_CUDA)
    # Create CUDA-enabled test
    add_executable(example_eikon_solver_cuda_test
        tests/example_eikon_solver_cuda_test.cpp
    )
    target_link_libraries(example_eikon_solver_cuda_test 
        eikon_solver_2d 
        hdf5_utils 
        Eigen3::Eigen 
        ${HDF5_LIBRARIES}
    )
endif()


# Create solver analysis executable


# Optional: Create a test executable (if you want to add unit tests)
# add_executable(test_eikon_solver test_eikon_solver.cpp)
# target_link_libraries(test_eikon_solver eikon_solver_2d Eigen3::Eigen)
