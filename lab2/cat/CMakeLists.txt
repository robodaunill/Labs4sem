cmake_minimum_required(VERSION 3.10)
project(CatEars)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libgmsh.so"  
    INTERFACE_INCLUDE_DIRECTORIES "/usr/include" 
)

add_executable(cat cat.cpp)
target_link_libraries(cat ${VTK_LIBRARIES} gmsh)