#include <iostream>
#include <set>
#include <gmsh.h>
#include <cmath>

int main(int argc, char **argv)
{
    gmsh::initialize();
    gmsh::model::add("cat");
    gmsh::merge("/home/robodaniil/Labs4sem/lab1/cat-v3.stl");

    std::vector<std::pair<int, int>> surfaceEntities;
    gmsh::model::getEntities(surfaceEntities, 2);

    std::vector<int> surfaceTags;
    for (auto entity : surfaceEntities) {
        surfaceTags.push_back(entity.second);
    }
    int surfaceLoopTag = gmsh::model::geo::addSurfaceLoop(surfaceTags);
    int volumeTag = gmsh::model::geo::addVolume({surfaceLoopTag});
    gmsh::model::geo::synchronize();

    double meshSize = 100.0;
    std::vector<std::pair<int, int>> pointEntities;
    gmsh::model::getEntities(pointEntities, 0);
    gmsh::model::mesh::setSize(pointEntities, meshSize);

    gmsh::option::setNumber("Mesh.MeshSizeFromPoints", 1);
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 1);
    
    gmsh::option::setNumber("Mesh.Algorithm3D", 1); 
    gmsh::option::setNumber("Mesh.Optimize", 1);
    gmsh::option::setNumber("Mesh.Smoothing", 4); 

    gmsh::model::mesh::generate(3);
    gmsh::fltk::run();
    gmsh::finalize();

    return 0;
}