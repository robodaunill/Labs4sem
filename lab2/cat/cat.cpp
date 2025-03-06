#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <tuple>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

using namespace std;

class CalcNode {
    friend class CalcMesh;

protected:
    double x;
    double y;
    double z;
    double smth;
    double vx;
    double vy;
    double vz;

public:
    CalcNode() : x(0.0), y(0.0), z(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0) {}

    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz)
        : x(x), y(y), z(z), smth(smth), vx(vx), vy(vy), vz(vz) {}

    void move(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }

    void setVelocity(double vx_, double vy_, double vz_) {
        vx = vx_;
        vy = vy_;
        vz = vz_;
    }

    std::tuple<double, double, double> getCoords() const {
        return std::make_tuple(x, y, z);
    }
};

class Element {
    friend class CalcMesh;
protected:
    unsigned long nodesIds[4];
};

class CalcMesh
{
protected:
    vector<CalcNode> nodes;
    vector<Element> elements;
    vector<unsigned long> tailNodes; // Узлы хвоста
    double tailCentroidX = 0.0;
    double tailCentroidZ = 0.0;

public:

    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) {
        nodes.resize(nodesCoords.size() / 3);
        for (unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
            double pointX = nodesCoords[i * 3];
            double pointY = nodesCoords[i * 3 + 1];
            double pointZ = nodesCoords[i * 3 + 2];
            double smth = pow(pointX, 2) + pow(pointY, 2) + pow(pointZ, 2); //Пример
            nodes[i] = CalcNode(pointX, pointY, pointZ, smth, 0.0, 0.0, 0.0);
        }

        elements.resize(tetrsPoints.size() / 4);
        for (unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
            elements[i].nodesIds[0] = tetrsPoints[i * 4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i * 4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i * 4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i * 4 + 3] - 1;
        }
    }

    void identifyTail() {
      tailCentroidX = 0.0;
      tailCentroidZ = 0.0;
        for (size_t i = 0; i < nodes.size(); ++i) {
            auto [x, y, z] = nodes[i].getCoords();
            if (y > 48.0) {
                tailNodes.push_back(i);
                tailCentroidX += x;
                tailCentroidZ += z;
            }
        }
        // Вычисляем центр хвоста
        if (!tailNodes.empty())
        {
              tailCentroidX /= tailNodes.size();
              tailCentroidZ /= tailNodes.size();
        }

        std::cout << "Found " << tailNodes.size() << " tail nodes." << std::endl;
        std::cout << "Tail centroid: (" << tailCentroidX << ", " << tailCentroidZ << ")" << std::endl;

    }


    void moveTail(double time) {
        const double omega = 2.0;  // Угловая скорость

        double adjustedTailCentroidZ = tailCentroidZ - 10.0;

        for (unsigned long nodeIndex : tailNodes) {
            auto [x, y, z] = nodes[nodeIndex].getCoords();

            // Вращение в плоскости XZ вокруг центра
            double dx = x - tailCentroidX;
            double dz = z - adjustedTailCentroidZ; // Используем adjustedTailCentroidZ
            double r = sqrt(dx * dx + dz * dz);
            double currentAngle = atan2(dz, dx);
            double newAngle = currentAngle + omega * time;
            double newX = tailCentroidX + r * cos(newAngle);
            double newZ = adjustedTailCentroidZ + r * sin(newAngle); // Используем adjustedTailCentroidZ

            double tau = 0.01;
            double vx = (newX - x) / tau;
            double vy = 0.0;
            double vz = (newZ - z) / tau;

            nodes[nodeIndex].setVelocity(vx, vy, vz);
        }
    }

    void doTimeStep(double tau) {
        for (unsigned int i = 0; i < nodes.size(); i++) {
            nodes[i].move(tau);
        }
    }

    void snapshot(unsigned int snap_number) {
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        for (unsigned int i = 0; i < nodes.size(); i++) {
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);
            double _vel[3] = { nodes[i].vx, nodes[i].vy, nodes[i].vz };
            vel->InsertNextTuple(_vel);
            smth->InsertNextValue(nodes[i].smth);
        }

        unstructuredGrid->SetPoints(dumpPoints);
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(smth);

        for (unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId(0, elements[i].nodesIds[0]);
            tetra->GetPointIds()->SetId(1, elements[i].nodesIds[1]);
            tetra->GetPointIds()->SetId(2, elements[i].nodesIds[2]);
            tetra->GetPointIds()->SetId(3, elements[i].nodesIds[3]);
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        string fileName = "cat-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main() {
    double tau = 0.01;
    int numSteps = 100;
    const unsigned int GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("cat");

    try {
        gmsh::merge("cat.stl");
    } catch(...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    double angle = 40;
    bool forceParametrizablePatches = true; 
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();

    std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto surf : s) sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});
    gmsh::model::geo::synchronize();

    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "4");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    gmsh::model::mesh::generate(3);

    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for (unsigned int i = 0; i < elementTypes.size(); i++) {
        if (elementTypes[i] != GMSH_TETR_CODE) continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if (tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " << nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    for (int i = 0; i < nodeTags.size(); ++i) {
        assert(i == nodeTags[i] - 1);
    }
    assert(tetrsNodesTags->size() % 4 == 0);

    CalcMesh mesh(nodesCoord, *tetrsNodesTags);

    // Идентификация и движение хвоста
    mesh.identifyTail();
    for (int step = 0; step < numSteps; ++step) {
        mesh.moveTail(step * tau);
        mesh.doTimeStep(tau);
        mesh.snapshot(step);
    }

    gmsh::finalize();
    return 0;
}