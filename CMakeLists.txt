cmake_minimum_required(VERSION 3.17)
project(SpatialCoherence)

set(CMAKE_CXX_STANDARD 20)

include_directories(/home/maarten/LAStools/LASzip/src /home/maarten/LAStools/LASlib/inc)

link_directories(/home/maarten/LAStools/LASlib/lib)

find_package(GDAL REQUIRED)

include_directories( include ${GDAL_INCLUDE_DIRS})

add_executable(SpatialCoherence main.cpp)

target_link_libraries(SpatialCoherence LASlib)
target_link_libraries(SpatialCoherence ${GDAL_LIBRARIES})

set_target_properties(SpatialCoherence PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
