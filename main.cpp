#include <iostream>
#include <stdexcept>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"

// #include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;


int main () {
	
	
		unsigned int p, q, b, c;
	
		// Input da terminale
		cout << "Inserisci il numero di lati per faccia (p): ";
		cin >> p;
	
		cout << "Inserisci il numero di facce per vertice (q): ";
		cin >> q;
	
		cout << "Inserisci il parametro di suddivisione (b): ";
		cin >> b;
	
		cout << "Inserisci la classe (c): ";
		cin >> c;
        
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