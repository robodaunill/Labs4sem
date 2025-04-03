#include "VTKWriter.hpp"
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkNew.h>

vtkSmartPointer<vtkPolyData> VTKWriter::createPolyData(const std::vector<Particle>& particles) {
    auto polydata = vtkSmartPointer<vtkPolyData>::New();
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto vertices = vtkSmartPointer<vtkCellArray>::New();

    // Добавляем точки
    for(const auto& p : particles) {
        vtkIdType id = points->InsertNextPoint(p.position.x, p.position.y, p.position.z);
        vertices->InsertNextCell(1, &id);
    }

    polydata->SetPoints(points);
    polydata->SetVerts(vertices);
    return polydata;
}

void VTKWriter::writeFrame(const std::vector<Particle>& particles,
                          const std::string& filename,
                          int frameNumber) {
    auto polydata = createPolyData(particles);
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    
    std::string fullPath = "../output/" + filename + "_" + 
                          std::to_string(frameNumber) + ".vtp";
    
    writer->SetFileName(fullPath.c_str());
    writer->SetInputData(polydata);
    writer->SetDataModeToBinary();
    
    if(writer->Write() == 0) {
        throw std::runtime_error("Failed to write VTK file: " + fullPath);
    }
}