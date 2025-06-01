#include <iostream>
#include <stdexcept>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"

#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;


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

		//i 4 file di output
        
        string outfilename0D="Cell0Ds.txt";
        string outfilename1D="Cell1Ds.txt";
        string outfilename2D="Cell2Ds.txt";
        string outfilename3D="Cell3Ds.txt";
		
		if (b==0 && c==0) throw runtime_error("Errore: uno tra b e c deve essere > 0.");

		//AGGIORNARE PER LA SECONDA CLASSE
		if(p==3 && q>=3 && q<=5 && (b == 0 || c == 0 || b == c)) {   //geodetico senza duale
			PolygonalMesh geodetico;
			if(b==0) {swap(b, c)};
			build_solido(p, q, b, c, geodetico);
			GeneraTuttiFile(geodetico, 
                     outfilename0D,
                     outfilename1D,
                     outfilename2D,
                     outfilename3D);
		} else if(q==3 && (p==4 || p==5) && (b == 0 || c == 0 || b == c)) {
			PolygonalMesh geodetico;
			PolygonalMesh duale;
			if(b==0) {swap(b, c)};
			build_solido(q, p, b, c, geodetico);
			build_duale(geodetico, duale);
			GeneraTuttiFile(duale, 
                     outfilename0D,
                     outfilename1D,
                     outfilename2D,
                     outfilename3D);
		} 
		else throw runtime_error("Valori NON accettabili!");

        
		/// Per visualizzare online le mesh:
		/// 1. Convertire i file .inp in file .vtu con https://meshconverter.it/it
		/// 2. Caricare il file .vtu su https://kitware.github.io/glance/app/

		Gedim::UCDUtilities utilities;
		{
			vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

			cell0Ds_properties[0].Label = "Marker";
			cell0Ds_properties[0].UnitLabel = "-";
			cell0Ds_properties[0].NumComponents = 1;
			
		/// DA TOGLIERE// 
			//vector<double> cell0Ds_marker(mesh.NumCell0Ds, 0.0);
			//for(const auto &m : mesh.MarkerCell0Ds)
			//	for(const unsigned int id: m.second)
			//		cell0Ds_marker.at(id) = m.first;

			//cell0Ds_properties[0].Data = cell0Ds_marker.data();

			utilities.ExportPoints("./Cell0Ds.inp",
								   mesh.Cell0DsCoordinates,
								   {});
		}

		{

			vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

			cell1Ds_properties[0].Label = "Marker";
			cell1Ds_properties[0].UnitLabel = "-";
			cell1Ds_properties[0].NumComponents = 1;
			
			/// DA TOGLIERE// 
			//vector<double> cell1Ds_marker(mesh.NumCell1Ds, 0.0);
			//for(const auto &m : mesh.MarkerCell1Ds)
			//	for(const unsigned int id: m.second)
			//		cell1Ds_marker.at(id) = m.first;

			//cell1Ds_properties[0].Data = cell1Ds_marker.data();

			utilities.ExportSegments("./Cell1Ds.inp",
									 mesh.Cell0DsCoordinates,
									 mesh.Cell1DsExtrema,
									 {},
									 {});
		}
		
	return 0;
	}