project(Rattleback CXX)
cmake_minimum_required(VERSION 2.8)

# Added the CMake files needed to find GSL and Eigen
set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Enable warnings
set(C_CPP_FLAGS "-Wall")

# Enable maximum optimizations for Release build
set(CMAKE_CXX_FLAGS_RELEASE "${C_CPP_FLAGS} -O3")

# Disable optimizations for debug build, enable maximum debugging symbols
set(CMAKE_CXX_FLAGS_DEBUG "${C_CPP_FLAGS} -O0 -ggdb3")

# Enable maximum optimizations for Release build, maximum debugging symbols
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${C_CPP_FLAGS} -O3 -ggdb3")

# Find required libraries
find_package(GSL 1.15 REQUIRED)
find_package(Eigen3 3.0 REQUIRED)
find_package(VTK REQUIRED)
add_definitions(-DVTK_LEGACY_REMOVE)
include(${VTK_USE_FILE})

include_directories(${EIGEN3_INCLUDE_DIR}
                    ${GSL_INCLUDE_DIR})

add_library(rattleback ellipsoid_no_slip.cc
                       paraboloid_no_slip.cc
                       ellipsoid_equilibrium.cc)

add_executable(ellipsoid_sim ellipsoid_sim.cc)
target_link_libraries(ellipsoid_sim
                      rattleback
                      ${GSL_LIBRARIES})

add_executable(paraboloid_sim paraboloid_sim.cc)
target_link_libraries(paraboloid_sim
                      rattleback
                      ${GSL_LIBRARIES})

add_executable(ellipsoid_equilibrium ellipsoid_eq.cc)
target_link_libraries(ellipsoid_equilibrium rattleback ${VTK_LIBRARIES})

# rattleback animation source files are missing:
#add_executable(RattleBackAnimation RattleBackAnimation.cc
#                                   RattleBackCallback.cc)
#target_link_libraries(RattleBackAnimation rattleback ${VTK_LIBRARIES})
