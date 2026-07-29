cmake_minimum_required(VERSION 3.18)
project(EikonSolver2D LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to enable CUDA support
option(ENABLE_CUDA "Enable CUDA GPU acceleration" OFF)
option(ENABLE_VTK "Enable VTK output support" OFF)
option(ENABLE_GMSH "Enable Gmsh output support" OFF)
option(BUILD_TESTS "Build unit tests" OFF)

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

# Optional: Find VTK for advanced output
if(ENABLE_VTK)
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
    add_definitions(-DVTK_ENABLED)
    message(STATUS "VTK support enabled")
endif()

# Add compiler flags for optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra")
if(ENABLE_CUDA)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 --use_fast_math --expt-relaxed-constexpr")
endif()

# Include directories
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(include)

# ============================================================
# SOURCE FILES
# ============================================================
set(EIKON_SOURCES
    src/eikon_solver_2d.cpp
)

# Add output utilities if they exist in separate files
if(EXISTS ${CMAKE_SOURCE_DIR}/src/vtk_output.cpp)
    list(APPEND EIKON_SOURCES src/vtk_output.cpp)
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/src/gmsh_output.cpp)
    list(APPEND EIKON_SOURCES src/gmsh_output.cpp)
endif()

if(ENABLE_CUDA)
    list(APPEND EIKON_SOURCES src/eikon_solver_cuda_impl.cu)
endif()

# ============================================================
# CREATE LIBRARIES
# ============================================================

# Main solver library
add_library(eikon_solver_2d ${EIKON_SOURCES})

# Set properties
if(ENABLE_CUDA)
    set_property(TARGET eikon_solver_2d PROPERTY CUDA_STANDARD 17)
    set_property(TARGET eikon_solver_2d PROPERTY CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES})
endif()

# Link dependencies
target_link_libraries(eikon_solver_2d 
    Eigen3::Eigen
)

if(ENABLE_CUDA)
    target_link_libraries(eikon_solver_2d CUDA::cudart)
endif()

if(ENABLE_VTK)
    target_link_libraries(eikon_solver_2d ${VTK_LIBRARIES})
endif()

# HDF5 utility library
add_library(hdf5_utils
    src/hdf5_utils.cpp
)

target_link_libraries(hdf5_utils 
    Eigen3::Eigen 
    ${HDF5_LIBRARIES}
)

# ============================================================
# EXAMPLES
# ============================================================

# CPU example
add_executable(example_eikon_solver 
    tests/example_eikon_solver.cpp
)
target_link_libraries(example_eikon_solver 
    eikon_solver_2d 
    hdf5_utils 
    Eigen3::Eigen 
    ${HDF5_LIBRARIES}
)
if(ENABLE_VTK)
    target_link_libraries(example_eikon_solver ${VTK_LIBRARIES})
endif()

# Triangle distance test (ADD THIS)
add_executable(triangle_distance 
    tests/triangle_distance.cpp
)
target_link_libraries(triangle_distance 
    eikon_solver_2d 
    hdf5_utils 
    Eigen3::Eigen 
    ${HDF5_LIBRARIES}
)
if(ENABLE_VTK)
    target_link_libraries(triangle_distance ${VTK_LIBRARIES})
endif()


# Triangle distance test (ADD THIS)
add_executable(test 
    tests/test.cpp
    src/polygon_utils.cpp      # ← ADD THESE THREE FILES
    src/cartesianMesh_utils.cpp         # ← ADD THESE THREE FILES
    src/io_utils.cpp           # ← ADD THESE THREE FILES
    src/prismMesh_utils.cpp
)
target_link_libraries(test
    eikon_solver_2d 
    hdf5_utils 
    Eigen3::Eigen 
    ${HDF5_LIBRARIES}
)
if(ENABLE_VTK)
    target_link_libraries(test ${VTK_LIBRARIES})
endif()


# CUDA example (if enabled)
if(ENABLE_CUDA)
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

# ============================================================
# COMPILER DEFINITIONS
# ============================================================
# Add definitions to the main library
target_compile_definitions(eikon_solver_2d PRIVATE 
    EIGEN_MPL2_ONLY
)

if(ENABLE_VTK)
    target_compile_definitions(eikon_solver_2d PRIVATE VTK_ENABLED)
endif()

if(ENABLE_GMSH)
    target_compile_definitions(eikon_solver_2d PRIVATE GMSH_ENABLED)
endif()

# ============================================================
# INSTALLATION
# ============================================================
install(TARGETS eikon_solver_2d hdf5_utils
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES 
    include/eikon_solver_2d.h
    include/hdf5_utils.h
    include/io_utils.hpp
    include/caresianMesh_utils.hpp
    include/polygon_utils.hpp
    include/prismMesh_tuils.hpp
    DESTINATION include
)

# ============================================================
# TESTING (optional)
# ============================================================
option(BUILD_TESTS "Build unit tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    
    add_executable(test_eikon_solver 
        tests/test_eikon_solver.cpp       
    )
    target_link_libraries(test_eikon_solver 
        eikon_solver_2d 
        GTest::GTest 
        GTest::Main
    )
    add_test(NAME EikonSolverTest COMMAND test_eikon_solver)
endif()

message(STATUS "Configuration complete!")


