#pragma once
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>
#include <vector>
#include "Simulation.hpp"

class VTKWriter {
public:
    static void writeFrame(const std::vector<Particle>& particles,
                          const std::string& filename,
                          int frameNumber);
private:
    static vtkSmartPointer<vtkPolyData> createPolyData(const std::vector<Particle>& particles);
};