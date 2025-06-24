#pragma once

#include <vector>
#include <map>
#include <list>
#include <Eigen/Dense>
// #include <queue>  TOGLIERE

/// Costruzione libreria PloyhedralLibrary
namespace PolyhedralLibrary {   
  struct PolyhedralMesh
  {
    /// Cell0D : vertici
    unsigned int Cell0D_num = 0;   // numero di Cell0D
    std::vector<unsigned int> Cell0D_id = {};   // vettore degli id di Cell0D : dim 1 x Cell0D_num
    Eigen::MatrixXd Cell0D_coordinate = {};   // matrice di coordinate di Cell0Dshape : X x Cell0D_num (x,y,z)
    std::map<unsigned int, std::list<unsigned int>> Cell0D_marker = {};   // mappa<marker, lista dei punti di quel marker>

    /// Cell1D : segmenti
    unsigned int Cell1D_num = 0;   // numero di Cell1D
    std::vector<unsigned int> Cell1D_id = {};   // vettore degli id di Cell1D : dim 1 x Cell1D_num
    Eigen::MatrixXi Cell1D_estremi = {};   // matrice degli id (in int) dei punti di estremi : dim 2 x Cell1D_num
    std::map<unsigned int, std::list<unsigned int>> Cell1D_marker = {};   // mappa<marker, lista dei punti di quel marker>

    /// Cell2D : facce
    unsigned int Cell2D_num = 0;   // numero di Cell2D
    std::vector<unsigned int> Cell2D_id = {};   // vettore degli id delle facce : dim 1 x Cell2D_num
    std::vector<std::vector<int>> Cell2D_vertici = {};   //  vettore di vettori degli id dei vertici di ciascuna faccia : dim 1 x Cell2D_num
    std::vector<int> Cell2D_numLati = {};   // vettore dei numeri di lati di ogni faccia : dim 1 x Cell2D_num
    std::vector<std::vector<int>> Cell2D_lati = {};   // vettore di vettori dei lati che compongono ciascuna faccia : dim 1 x CEll2D_num
    
    /// Cell3D : solidi
    unsigned int Cell3D_num = 0;   // numero di Cell3D
    std::vector<unsigned int> Cell3D_id = {};   // vettore degli id dei solidi : dim 1 x Cell3D_num
    std::vector<std::vector<int>> Cell3D_vertici = {};   //  vettore di vettori degli id dei vertici di ciascuna solido : dim 1 x Cell3D_num
    std::vector<std::vector<int>> Cell3D_lati = {};   // vettore di vettori dei lati che compongono ciascun solido : dim 1 x Cell3D_num
    std::vector<std::vector<int>> Cell3D_facce = {};   // vettore di vettori di facce che compongono ciascun solido : dim 1 x Cell3D_num
    
	/* === Cammino Minimo Path === TOGLIERE
  unsigned int num_archiPath = 0;
  double lunghezza_Path = 0.0;
  */
  };
} 