#always start with the minimum required version
cmake_minimum_required(VERSION 3.16)

#set the project name
project(CajeteProto VERSION 0.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#set the C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17) 

add_compile_options(-O3)

# Add the serial host mode by default
# Currently turning it of would break the code
option(ENABLE_SERIAL_HOST "Enable Serial Host Mode" ON)
if(ENABLE_SERIAL_HOST)
	add_definitions(-DCAJETE_ENABLE_HOST_SERIAL)
endif()

# Add debugging option
option(USE_DEBUG "Enter debug mode" ON)
if (USE_DEBUG)
	add_definitions(-DCAJETE_DEBUG)
endif()

# In the build script simpy do this:
#
# Set the variable SUNDIALS_DIR to the SUNDIALS instdir.
# When using the cmake CLI command, this can be done like so:
#   cmake -D SUNDIALS_DIR=/path/to/sundials/installation

#find the sundials package for ODE solving 
find_package(SUNDIALS REQUIRED)

if(SUNDIALS_FOUND)
	option(SUNDIALS_FOUND ON)
endif() 

include_directories(include)
include_directories(thirdparty/simdjson)
include_directories(src)
include_directories(src/PlantModel)
include_directories(src/Grids)
include_directories(src/CellComplexes)
include_directories(src/CellComplexes/Utils)
include_directories(src/Solvers)
add_subdirectory(thirdparty/vtu11)
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(tools)
add_subdirectory(simulators)

set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/simulators/configuration")
set(DESTINATION_DIR "${CMAKE_BINARY_DIR}/simulators")
file(COPY ${SOURCE_DIR} DESTINATION ${DESTINATION_DIR})
