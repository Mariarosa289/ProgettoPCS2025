#include <iostream>
#include <stdexcept>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"

// #include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;


int main () {
		unsigned int p = 3;   // Numero lati per faccia (solo 3 supportato)
        unsigned int q = 5;   // Numero facce per vertice (5 = icosaedro)
        unsigned int b = 1;   // Parametro di suddivisione
        unsigned int c = 0;   // Solo classe I (c=0)
        
        PolygonalMesh mesh;
		
        // Costruzione della mesh
        build_solido(p, q, b, c, mesh);
        
        //i 4 file di output
        
        string outfilename0D="Cell0Ds.txt";
        string outfilename1D="Cell1Ds.txt";
        string outfilename2D="Cell2Ds.txt";
        string outfilename3D="Cell3Ds.txt";
        
        // generiamo i 4 file di output
        
		GeneraTuttiFile(mesh, 
                     outfilename0D,
                     outfilename1D,
                     outfilename2D,
                     outfilename3D);		
		
	return 0;
	}