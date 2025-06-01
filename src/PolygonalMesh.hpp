#pragma once

#include <vector>
#include <map>
#include <list>
#include <Eigen/Eigen>

namespace PolygonalLibrary {

struct PolygonalMesh
{
    // === Cell0D ===
    unsigned int NumCell0Ds = 0; ///< number of Cell0D
    std::vector<unsigned int> Cell0DsId = {}; ///< Cell0D ids, size NumCell0Ds
    Eigen::MatrixXd Cell0DsCoordinates = {}; ///< Cell0D coordinates, shape: X x NumCell0Ds (x,y)
    std::map<unsigned int, std::list<unsigned int>> MarkerCell0Ds = {}; ///< Cell0D markers, e.g., ShortestPath marker

    // === Cell1D ===
    unsigned int NumCell1Ds = 0; ///< number of Cell1D
    std::vector<unsigned int> Cell1DsId = {}; ///< Cell1D ids, size NumCell1Ds
    Eigen::MatrixXi Cell1DsExtrema = {}; ///< Cell1D extrema (from, to), shape: X x NumCell1Ds
    std::map<unsigned int, std::list<unsigned int>> MarkerCell1Ds = {}; ///< Cell1D markers, e.g., ShortestPath marker

    // === Cell2D ===
    unsigned int NumCell2Ds = 0; ///< number of Cell2D
    std::vector<unsigned int> Cell2DsId = {}; ///< Cell2D ids, size NumCell2Ds
    std::vector<std::vector<int>> Cell2DsVertices = {}; ///< Cell2D vertex ids, variable size per cell
    std::vector<int> Cell2DsNumEdges = {}; ///< Number of edges per Cell2D, size NumCell2Ds
    std::vector<std::vector<int>> Cell2DsEdges = {}; ///< Cell2D edge ids (Cell1D), variable size per cell
    std::map<unsigned int, std::list<unsigned int>> MarkerCell2Ds = {}; ///< Cell2D markers

  
  // === Cell3D ===
    unsigned int NumCell3Ds = 0;
    std::vector<unsigned int> Cell3DsId = {};
    std::vector<std::vector<int>> Cell3DsVertices = {};
    std::vector<std::vector<int>> Cell3DsEdges = {};
    std::vector<std::vector<int>> Cell3DsFaces = {};
    std::map<unsigned int, std::list<unsigned int>> MarkerCell3Ds = {};
	
};

} // namespace PolygonalLibrary
