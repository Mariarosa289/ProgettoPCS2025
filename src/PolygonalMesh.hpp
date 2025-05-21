#pragma once
#include <vector>
#include <array>

struct Vertex {
    unsigned int id;
    std::array<double, 3> coordinates;
    bool ShortPath = false;
};

struct Edge {
    unsigned int id;
    unsigned int origin, end;
    bool ShortPath = false;
};

struct Face {
    unsigned int id;
    std::vector<unsigned int> vertices;
    std::vector<unsigned int> edges;
};

struct Polyhedron {
    unsigned int id;
    std::vector<unsigned int> vertices;
    std::vector<unsigned int> edges;
    std::vector<unsigned int> faces;
};

struct PolygonalMesh {
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::vector<Polyhedron> polyhedra;
};

