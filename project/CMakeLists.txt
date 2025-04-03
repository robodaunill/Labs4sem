cmake_minimum_required(VERSION 3.10)
project(DropletSimulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Vulkan REQUIRED)
find_package(VTK REQUIRED COMPONENTS
  CommonCore
  CommonDataModel
  IOXML
)

add_executable(droplet_sim
    src/main.cpp
    src/Simulation.cpp
    src/VTKWriter.cpp
)

target_include_directories(droplet_sim PRIVATE 
    ${VTK_INCLUDE_DIRS}
    include
)

target_link_libraries(droplet_sim
    Vulkan::Vulkan
    ${VTK_LIBRARIES}
)