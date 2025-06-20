#pragma once
#include "PolyhedralMesh.hpp"
#include <string>
#include <map>
#include <vector>
#include <Eigen/Eigen>
#include <queue>
#include <list>

using namespace std;
using namespace PolyhedralLibrary;
using namespace Eigen;




void build_classe_1(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolyhedralMesh& mesh);

void build_classe_2(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolyhedralMesh& mesh);

pair<int, bool> get_or_add_edge(
    unsigned int a, unsigned int b,
    map<pair<unsigned int, unsigned int>, unsigned int>& map1D,
    PolyhedralMesh& mesh);

array<double,3> PuntoMedio(const array<double,3>& A, const array<double,3>& B);

array<double,3> Baricentro(const array<double,3>& A, const array<double,3>& B, const array<double,3>& C);

unsigned int salva_vertice_norm(const array<double,3>& coord,
                           map<array<double,3>, unsigned int>& map0D,
                           PolyhedralMesh& mesh,
                           unsigned int& vid) ;

void build_tetra(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce);

void build_octa(vector<array<double, 3>>& vertici,vector<array<int, 3>>& facce);

void build_ico(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce);

array<double, 3> Normalizza(const array<double, 3>& p);

bool GeneraFileCell0D(const PolyhedralMesh& mesh, const string& outfilename);

bool GeneraFileCell1D(const PolyhedralMesh& mesh, const string& outfilename);

bool GeneraFileCell2D(const PolyhedralMesh& mesh, const string& outfilename);

bool GeneraFileCell3D(const PolyhedralMesh& mesh, const string& outfilename);

void GeneraTuttiFile(const PolyhedralMesh& mesh, 
                     const string& outfilename0D,
                     const string& outfilename1D,
                     const string& outfilename2D,
                     const string& outfilename3D);

static Vector3d CentroideCoordinate(const vector<int>& VerticesId, const MatrixXd& coordinates);

void build_duale(const PolyhedralMesh& geodetico, PolyhedralMesh& Goldberg);

void Dijkstra(PolyhedralMesh& mesh, unsigned int vertice_iniziale, unsigned int vertice_finale);




