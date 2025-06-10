#pragma once
#include "PolygonalMesh.hpp"
#include <string>
#include <map>
#include <vector>
#include <Eigen/Eigen>
#include <queue>
#include <list>

using namespace std;
using namespace PolygonalLibrary;
using namespace Eigen;




void build_solido(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolygonalMesh& mesh);

unsigned int salva_vertice_norm(const array<double,3>& coord,
                           map<array<double,3>, unsigned int>& map0D,
                           PolygonalMesh& mesh,
                           unsigned int& vid) ;

void build_tetra(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce);

void build_octa(vector<array<double, 3>>& vertici,vector<array<int, 3>>& facce);

void build_ico(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce);

bool ControllaInput(unsigned int p, unsigned int q, unsigned int b, unsigned int c);

array<double, 3> Normalizza(const array<double, 3>& p);

bool GeneraFileCell0D(const PolygonalMesh& mesh, const string& outfilename);

bool GeneraFileCell1D(const PolygonalMesh& mesh, const string& outfilename);

bool GeneraFileCell2D(const PolygonalMesh& mesh, const string& outfilename);

bool GeneraFileCell3D(const PolygonalMesh& mesh, const string& outfilename);

void GeneraTuttiFile(const PolygonalMesh& mesh, 
                     const string& outfilename0D,
                     const string& outfilename1D,
                     const string& outfilename2D,
                     const string& outfilename3D);

static Vector3d CentroideCoordinate(const vector<int>& VerticesId, const MatrixXd& coordinates);

void build_duale(const PolygonalMesh& geodetico, PolygonalMesh& Goldberg);

void Dijkstra(PolygonalMesh& mesh, unsigned int vertice_iniziale, unsigned int vertice_finale);

void triangolazione_baricentrica_b_c(
    const array<double,3>& A,
    const array<double,3>& B,
    const array<double,3>& C,
    unsigned int b,
    map<array<double,3>, unsigned int>& map0D,
    map<pair<unsigned int,unsigned int>, unsigned int>& map1D,
    PolygonalMesh& mesh,
    unsigned int& fid
) ;


