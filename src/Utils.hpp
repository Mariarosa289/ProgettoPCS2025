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


array<double, 3> normalizza(const array<double, 3>& p);

array<double, 3> arrotonda_punto(const array<double, 3>& p);

void build_tetra(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce);

void build_octa(vector<array<double, 3>>& vertici,vector<array<int, 3>>& facce);

void build_ico(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce);

unsigned int salva_vertice_norm(const array<double,3>& coord,
                                map<array<double,3>, unsigned int>& map0D,
                                PolyhedralMesh& mesh,
                                unsigned int& vid) ;

array<double,3> calcola_baricentro(const array<double,3>& A, 
                                   const array<double,3>& B, 
                                   const array<double,3>& C);

array<double,3> calcola_puntomedio(const array<double,3>& A, 
                                   const array<double,3>& B); 

void build_classe_1(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolyhedralMesh& mesh);

void build_classe_2(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolyhedralMesh& mesh);

bool genera_file_Cell0D(const PolyhedralMesh& mesh, const string& outfilename);

bool genera_file_Cell1D(const PolyhedralMesh& mesh, const string& outfilename);

bool genera_file_Cell2D(const PolyhedralMesh& mesh, const string& outfilename);

bool genera_file_Cell3D(const PolyhedralMesh& mesh, const string& outfilename);

void genera_tutti_file(const PolyhedralMesh& mesh, 
                       const string& outfilename0D,
                       const string& outfilename1D,
                       const string& outfilename2D,
                       const string& outfilename3D);

static Vector3d calcola_coordinate_centroide(const vector<int>& VerticesId, const MatrixXd& Cell0D_coordinate);

void build_duale(const PolyhedralMesh& geodetico, PolyhedralMesh& Goldberg);

void Dijkstra(PolyhedralMesh& mesh, unsigned int vertice_iniziale, unsigned int vertice_finale);




