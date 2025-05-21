#pragma once
#include "PolygonalMesh.hpp"
#include <string>

using namespace std;



void GenerateGeodesicPolyhedron(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolygonalMesh &mesh);

void GenerateDual(const PolygonalMesh &original, PolygonalMesh &dual);

void WriteMeshToTxt(const PolygonalMesh &mesh, const std::string &prefix);

void ProjectVerticesOnSphere(PolygonalMesh &mesh);

bool ComputeShortestPath(PolygonalMesh &mesh, unsigned int startId, unsigned int endId, std::vector<unsigned int> &path, double &totalLength);

bool CheckInput(unsigned int p, unsigned int q, unsigned int b, unsigned int c);

//**************************************************************************** */

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

