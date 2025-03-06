cmake_minimum_required(VERSION 3.10)
project(Tetr3D)

add_executable(tetr3d tetr3d.cpp)

add_library(gmsh SHARED IMPORTED)  # Создаем imported target для gmsh
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libgmsh.so"  
    INTERFACE_INCLUDE_DIRECTORIES "/usr/include" 
)


# Прилинковываем VTK и GMSH
target_link_libraries(tetr3d ${VTK_LIBRARIES} gmsh) # Используем имя target-а gmsh