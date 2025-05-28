#pragma once
#include "PolygonalMesh.hpp"
#include <string>
#include <map>

using namespace std;




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

/// SCRIVERE CODICI PER QUESTE TRE FUNZIONI
void GenerateDual(const PolygonalMesh &original, PolygonalMesh &dual);

void WriteMeshToTxt(const PolygonalMesh &mesh, const std::string &prefix);

bool ComputeShortestPath(PolygonalMesh &mesh, unsigned int startId, unsigned int endId, std::vector<unsigned int> &path, double &totalLength);



//**************************************************************************** */

/// RIVEDERE COSA SONO

namespace PolygonalLibrary
{
/// Import the triangular mesh and test if the mesh is correct
/// mesh: a TriangularMesh struct
/// return the result of the reading, true if is success, false otherwise
bool ImportMesh(PolygonalMesh& mesh);   // RIVEDERE

/// Import the Cell0D properties from Cell0Ds.csv file
/// mesh: a TriangularMesh struct
/// return the result of the reading, true if is success, false otherwise
bool build_tetra(PolygonalMesh& mesh);

/// Import the Cell1D properties from Cell1Ds.csv file
/// mesh: a TriangularMesh struct
/// return the result of the reading, true if is success, false otherwise
bool build_octa(PolygonalMesh& mesh);

/// Import the Cell2D properties from Cell2Ds.csv file
/// mesh: a TriangularMesh struct
/// return the result of the reading, true if is success, false otherwise
bool build_ico(PolygonalMesh& mesh);


}

