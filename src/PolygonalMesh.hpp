// RIVEDERE LE #include

#pragma once
#include <vector>
#include <array>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

struct Cell0D {
    unsigned int id;
    array<double, 3> coordinate;
    bool ShortPath = false;   // SCRIVI CODICE DEL CAMMINO MINIMO 
};

struct Cell1D {
    unsigned int id;
    unsigned int origine, fine;
    bool ShortPath = false;   // SCRIVI CODICE DEL CAMMINO MINIMO
};

struct Cell2D {
    unsigned int id;
    vector<unsigned int> vertici;
    vector<unsigned int> spigoli;
};

struct Cell3D {
    unsigned int id;
    vector<unsigned int> vertici;
    vector<unsigned int> spigoli;
    vector<unsigned int> facce;
};

struct PolygonalMesh {
    vector<Cell0D> vertici;
    vector<Cell1D> spigoli;
    vector<Cell2D> facce;
    vector<Cell3D> poliedri;
};

