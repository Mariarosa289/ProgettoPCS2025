// RIVEDERE LE #include

#pragma once
#include <vector>
#include <array>

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
    vector<unsigned int> Cell0D;
    vector<unsigned int> Cell1D;
};

struct Cell3D {
    unsigned int id;
    vector<unsigned int> Cell0D;
    vector<unsigned int> Cell1D;
    vector<unsigned int> Cell2D;
};

struct PolygonalMesh {
    vector<Cell0D> Cell0D;
    vector<Cell1D> Cell1D;
    vector<Cell2D> Cell2D;
    vector<Cell3D> Cell3D;
};

