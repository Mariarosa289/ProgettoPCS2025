#include <iostream>
#include <stdexcept>
#include <vector>
#include <Eigen/Dense>
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include <queue>
#include <map>
#include <list>

#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;


int main () {
	
	
		unsigned int p, q, b, c, vertice_iniziale, vertice_finale;
		auto r=0;
	
		// Input da terminale
		cout << "Inserisci il numero di lati per faccia (p): ";
		cin >> p;
	
		cout << "Inserisci il numero di facce per vertice (q): ";
		cin >> q;
	
		cout << "Inserisci il parametro di suddivisione (b): ";
		cin >> b;
	
		cout << "Inserisci la classe (c): ";
		cin >> c;

		if (b==0 && c==0) throw runtime_error("Errore: uno tra b e c deve essere > 0.");

		PolygonalMesh mesh;
		
		if(p==3 && q>=3 && q<=5 && (b == 0 || c == 0)) {   //geodetico in classe_1 senza duale	
			if(b==0) {swap(b, c);}
			build_classe_1(p, q, b, c, mesh);	

		} else if(q==3 && (p==4 || p==5) && (b == 0 || c == 0)) {   //duale in classe_1
			PolygonalMesh geodetico;

			if(b==0) {swap(b, c);}
			build_classe_1(q, p, b, c, geodetico);
			build_duale(geodetico, mesh);
	
		} else if(p==3 && q>=3 && q<=5 && (b==c)) {   //geodetico in classe_2 senza duale
			build_classe_2(p, q, b, c, mesh);

		} else ifelse if(q==3 && (p==4 || p==5) && (b == c)) {   // duale in classe_2
			PolygonalMesh geodetico;

			build_classe_2(q, p, b, c, geodetico);
			build_duale(geodetico, mesh);

		} else throw runtime_error("Valori NON accettabili!");


		
		/// chiede per il cammino minimo
		cout << "Digita 1 se vuoi calcolare il cammino minimo, altrimenti digita un qualsiasi carattere: ";
		cin >> r;

		if (r==1) {
			cout << "Inserisci l'ID del vertice iniziale: ";
			cin >> vertice_iniziale;
			cout << "Inserisci l'ID del vertice finale: ";
			cin >> vertice_finale;

			if (vertice_finale <0 || vertice_finale >mesh.NumCell0Ds-1 || vertice_iniziale <0 || vertice_iniziale >mesh.NumCell0Ds-1  ) 
			throw runtime_error("ID dei vertici inseriti non validi. Inserire di nuovo tutti i parametri.");

			Dijkstra(mesh, vertice_iniziale, vertice_finale);
			cout << "Il cammino minimo che va da " << vertice_iniziale << " a " << vertice_finale << " è formato da " 
			<< mesh.num_archiPath << " lati e ha lunghezza " << mesh.lunghezza_Path << "." << endl;
		}

		//i 4 file di output
        
        string outfilename0D="Cell0Ds.txt";
        string outfilename1D="Cell1Ds.txt";
        string outfilename2D="Cell2Ds.txt";
        string outfilename3D="Cell3Ds.txt";
		
		std::cout << "NumCell0Ds = " << mesh.NumCell0Ds << ", size Cell0DsId = " << mesh.Cell0DsId.size() << ", Cell0DsCoordinates.cols() = " << mesh.Cell0DsCoordinates.cols() << std::endl;
		std::cout << "NumCell1Ds = " << mesh.NumCell1Ds << ", size Cell1DsId = " << mesh.Cell1DsId.size() << ", Cell1DsExtrema.cols() = " << mesh.Cell1DsExtrema.cols() << std::endl;
		std::cout << "NumCell2Ds = " << mesh.NumCell2Ds << ", size Cell2DsId = " << mesh.Cell2DsId.size() << std::endl;
		
		std::cout << "Cell2DsId.size() = " << mesh.Cell2DsId.size() << std::endl;
		std::cout << "Cell2DsVertices.size() = " << mesh.Cell2DsVertices.size() << std::endl;
		std::cout << "Cell2DsNumEdges.size() = " << mesh.Cell2DsNumEdges.size() << std::endl;
		
		GeneraTuttiFile(mesh, 
                     outfilename0D,
                     outfilename1D,
                     outfilename2D,
                     outfilename3D);

        
		/// Per visualizzare online le mesh:
		/// 1. Convertire i file .inp in file .vtu con https://meshconverter.it/it
		/// 2. Caricare il file .vtu su https://kitware.github.io/glance/app/

		Gedim::UCDUtilities utilities;
		{
			vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

			cell0Ds_properties[0].Label = "Marker";
			cell0Ds_properties[0].UnitLabel = "-";
			cell0Ds_properties[0].NumComponents = 1;
			
		
			vector<double> cell0Ds_marker(mesh.NumCell0Ds, 0.0);
			for(const auto& m : mesh.MarkerCell0Ds)
				for(const unsigned int id: m.second)
					cell0Ds_marker.at(id) = m.first;

			cell0Ds_properties[0].Data = cell0Ds_marker.data();

			utilities.ExportPoints("./Cell0Ds.inp",
								   mesh.Cell0DsCoordinates,
								   cell0Ds_properties);
		}

		{

			vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

			cell1Ds_properties[0].Label = "Marker";
			cell1Ds_properties[0].UnitLabel = "-";
			cell1Ds_properties[0].NumComponents = 1;
			
			
			vector<double> cell1Ds_marker(mesh.NumCell1Ds, 0.0);
			for(const auto &m : mesh.MarkerCell1Ds)
				for(const unsigned int id: m.second)
					cell1Ds_marker.at(id) = m.first;

			cell1Ds_properties[0].Data = cell1Ds_marker.data();

			utilities.ExportSegments("./Cell1Ds.inp",
									 mesh.Cell0DsCoordinates,
									 mesh.Cell1DsExtrema,
									 {},
									 cell1Ds_properties);
		}
		
	return 0;
	}
	
	

