#include <iostream>
#include <stdexcept>
#include <vector>
#include <Eigen/Eigen>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include <queue>
#include <map>
#include <list>

#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;


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
		
		//AGGIORNARE PER LA SECONDA CLASSE
		if(p==3 && q>=3 && q<=5 && (b == 0 || c == 0 || b == c)) {   //geodetico senza duale

			if(b==0) {swap(b, c);}
			build_solido(p, q, b, c, mesh);	
		} else if(q==3 && (p==4 || p==5) && (b == 0 || c == 0 || b == c)) {
			PolygonalMesh geodetico;
			
			if(b==0) {swap(b, c);}
			build_solido(q, p, b, c, geodetico);
			build_duale(geodetico, mesh);
	
		} 
		else throw runtime_error("Valori NON accettabili!");
		
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