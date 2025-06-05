#include "Utils.hpp"
#include "PolygonalMesh.hpp"
#include <vector>
#include <Eigen/Eigen>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <tuple>
#include <iomanip>
#include <cmath>
#include <numeric> 
#include <queue>
#include <list>

using namespace PolygonalLibrary;
using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------*/

/*
Funzione per normalizzare le coordinate di un punto, ovvero la proiezione dei punti sulla sfera
- function: Normalizza
- input 1: array con le 3 coordinate non normalizzate del punto p
- return: array con le 3 coordinate normalizzate del punto p
*/

array<double, 3> Normalizza(const array<double, 3>& p) {
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);   // calcolo la norma per normalizzare
    return {p[0]/norm, p[1]/norm, p[2]/norm};   
}

/*-----------------------------------------------------------------------------------------------*/

array<double, 3> round_point(const array<double, 3>& p, double eps = 1e-6) {
    return {
        round(p[0] / eps) * eps,
        round(p[1] / eps) * eps,
        round(p[2] / eps) * eps
    };
}

/*-----------------------------------------------------------------------------------------------*/

bool ControllaInput(unsigned int p, unsigned int q, unsigned int b, unsigned int c) {
    return p >= 3 && q >= 3 && q <= 5 && (b == 0 || c == 0 || b == c);
}

/*-----------------------------------------------------------------------------------------------*/

void build_ico(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce) {
    const double phi = (1.0 + sqrt(5.0)) / 2.0;
    vertici = {
        {-1,  phi,  0}, { 1,  phi,  0}, {-1, -phi,  0}, { 1, -phi,  0},
        { 0, -1,  phi}, { 0,  1,  phi}, { 0, -1, -phi}, { 0,  1, -phi},
        { phi,  0, -1}, { phi,  0,  1}, {-phi,  0, -1}, {-phi,  0,  1}
    };
    facce = {
        {0,11,5}, {0,5,1}, {0,1,7}, {0,7,10}, {0,10,11},
        {1,5,9}, {5,11,4}, {11,10,2}, {10,7,6}, {7,1,8},
        {3,9,4}, {3,4,2}, {3,2,6}, {3,6,8}, {3,8,9},
        {4,9,5}, {2,4,11}, {6,2,10}, {8,6,7}, {9,8,1}
    };
}

/*-----------------------------------------------------------------------------------------------*/

void build_octa(vector<array<double, 3>>& vertici,vector<array<int, 3>>& facce) {
    vertici = {
        { 1,  0,  0}, {-1,  0,  0}, { 0,  1,  0},
        { 0, -1,  0}, { 0,  0,  1}, { 0,  0, -1}
    };
    facce = {
        {0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
        {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}
    };
}

/*-----------------------------------------------------------------------------------------------*/

void build_tetra(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce) {
    vertici = {
        {  1,  1,  1 },
        { -1, -1,  1 },
        { -1,  1, -1 },
        {  1, -1, -1 }
    };
    facce = {
        {0, 1, 2}, {0, 3, 1}, {0, 2, 3}, {1, 3, 2}
    };
}

/*-----------------------------------------------------------------------------------------------*/

unsigned int salva_vertice_norm(const array<double,3>& coord,
                           map<array<double,3>, unsigned int>& map0D,
                           PolygonalMesh& mesh) 
{
    static unsigned int vid = 0;
    auto norm = Normalizza(coord);
    auto rounded = round_point(norm);
    auto it = map0D.find(rounded);
    if (it != map0D.end()) return it->second;
    map0D[rounded] = vid;
    mesh.Cell0DsId.push_back(vid);
    mesh.Cell0DsCoordinates.conservativeResize(3, vid + 1);
    for (int i = 0; i < 3; ++i)
        mesh.Cell0DsCoordinates(i, vid) = norm[i];
    return vid++;
}

/*-----------------------------------------------------------------------------------------------*/

void build_solido(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolygonalMesh& mesh) {
	
	if (!ControllaInput(p, q, b, c)) throw runtime_error("Input non valido");   
    if (p != 3 || c != 0) throw runtime_error("Supportati solo solidi di classe I (p=3, c=0)");   //MAIN:swap b e c
    
    vector<array<double, 3>> vertici_iniziali;
    vector<array<int, 3>> facce_iniziali;

    if (q == 5) build_ico(vertici_iniziali, facce_iniziali);
    else if (q == 4) build_octa(vertici_iniziali, facce_iniziali);
    else build_tetra(vertici_iniziali, facce_iniziali);

    map<array<double, 3>, unsigned int> map0D;
    map<pair<unsigned int,unsigned int>, unsigned int> map1D;
    unsigned int fid = 0;

    for (const auto& face : facce_iniziali) {
        const auto& A = vertici_iniziali[face[0]];
        const auto& B = vertici_iniziali[face[1]];
        const auto& C = vertici_iniziali[face[2]];

        vector<vector<unsigned int>> griglia(b + 1);
        for (unsigned int i = 0; i <= b; ++i) {
            griglia[i].resize(i + 1);
            for (unsigned int j = 0; j <= i; ++j) {
                double u = 1.0 - static_cast<double>(i)/b;
                double v = static_cast<double>(i - j)/b;
                double w = static_cast<double>(j)/b;
                array<double,3> P = {
                    u*A[0] + v*B[0] + w*C[0],
                    u*A[1] + v*B[1] + w*C[1],
                    u*A[2] + v*B[2] + w*C[2]
                };
                griglia[i][j] = salva_vertice_norm(P, map0D, mesh);
            }
        }

        for (unsigned int i = 0; i < b; ++i) {
            for (unsigned int j = 0; j < i + 1; ++j) {
                unsigned int v1 = griglia[i][j];
                unsigned int v2 = griglia[i+1][j];
                unsigned int v3 = griglia[i+1][j+1];
                mesh.Cell2DsId.push_back(fid++);
                mesh.Cell2DsVertices.push_back({(int)v1, (int)v2, (int)v3});
                vector<int> edge_ids;
                for (int k = 0; k < 3; ++k) {
                    unsigned int a = mesh.Cell2DsVertices.back()[k];
                    unsigned int b = mesh.Cell2DsVertices.back()[(k+1)%3];
                    if (a > b) swap(a, b);
                    auto key = make_pair(a, b);
                    if (!map1D.count(key)) {
                        map1D[key] = mesh.Cell1DsId.size();
                        mesh.Cell1DsId.push_back(map1D[key]);
                        mesh.Cell1DsExtrema.conservativeResize(2, map1D[key] + 1);
                        mesh.Cell1DsExtrema(0, map1D[key]) = a;
                        mesh.Cell1DsExtrema(1, map1D[key]) = b;
                    }
                    edge_ids.push_back(map1D[key]);
                }
                mesh.Cell2DsEdges.push_back(edge_ids);
                mesh.Cell2DsNumEdges.push_back(3);

                if (j < i) {
                    unsigned int v4 = griglia[i][j+1];
                    mesh.Cell2DsId.push_back(fid++);
                    mesh.Cell2DsVertices.push_back({(int)v1, (int)v3, (int)v4});
                    vector<int> edge_ids2;
                    for (int k = 0; k < 3; ++k) {
                        unsigned int a = mesh.Cell2DsVertices.back()[k];
                        unsigned int b = mesh.Cell2DsVertices.back()[(k+1)%3];
                        if (a > b) swap(a, b);
                        auto key = make_pair(a, b);
                        if (!map1D.count(key)) {
                            map1D[key] = mesh.Cell1DsId.size();
                            mesh.Cell1DsId.push_back(map1D[key]);
                            mesh.Cell1DsExtrema.conservativeResize(2, map1D[key] + 1);
                            mesh.Cell1DsExtrema(0, map1D[key]) = a;
                            mesh.Cell1DsExtrema(1, map1D[key]) = b;
                        }
                        edge_ids2.push_back(map1D[key]);
                    }
                    mesh.Cell2DsEdges.push_back(edge_ids2);
                    mesh.Cell2DsNumEdges.push_back(3);
                }
            }
        }
    }

    mesh.NumCell0Ds = mesh.Cell0DsId.size();
    mesh.NumCell1Ds = mesh.Cell1DsId.size();
    mesh.NumCell2Ds = mesh.Cell2DsId.size();


	// Salva la Cell3D (unica)
    mesh.NumCell3Ds = 1;
    mesh.Cell3DsId.push_back(0);

    // Tutti i vertici usati nella mesh
    vector<int> cell3D_vertices(mesh.NumCell0Ds);
    iota(cell3D_vertices.begin(), cell3D_vertices.end(), 0);
    mesh.Cell3DsVertices.push_back(cell3D_vertices);

    // Tutti gli spigoli usati nella mesh
    vector<int> cell3D_edges(mesh.NumCell1Ds);
    iota(cell3D_edges.begin(), cell3D_edges.end(), 0);
    mesh.Cell3DsEdges.push_back(cell3D_edges);

    // Tutte le facce usate nella mesh
    vector<int> cell3D_faces(mesh.NumCell2Ds);
    iota(cell3D_faces.begin(), cell3D_faces.end(), 0);
    mesh.Cell3DsFaces.push_back(cell3D_faces);


	
    mesh.Cell2DsId.shrink_to_fit(); // shrink_to_fit libera la memoria inutilizzata
    mesh.Cell2DsEdges.shrink_to_fit();
    mesh.Cell2DsVertices.shrink_to_fit();
    mesh.Cell2DsNumEdges.shrink_to_fit();

    
}

/*-----------------------------------------------------------------------------------------------*/

/// Funzione che genera tutti i file
void GeneraTuttiFile(const PolygonalMesh& mesh, 
                     const string& outfilename0D,
                     const string& outfilename1D,
                     const string& outfilename2D,
                     const string& outfilename3D) {
	 if (! GeneraFileCell0D (mesh, outfilename0D))
	 {
		 cerr<<"Problemi di esportazione per Cell0Ds.txt"<< endl;
		 
	  }
	  else 
	  	cout<<"Esportazione terminata con successo per Cell0Ds.txt"<< endl;
	  	
	  
	  if (! GeneraFileCell1D (mesh, outfilename1D))
	 {
		 cerr<<"Problemi di esportazione per Cell1Ds.txt"<< endl;
		 
	  }
	  else
	  	cout<<"Esportazione terminata con successo per Cell1Ds.txt"<< endl;
	  	
	  	
	  if (! GeneraFileCell2D (mesh, outfilename2D))
	 {
		 cerr<<"Problemi di esportazione per Cell2Ds.txt"<< endl;
		 
	  }
	  else
	  	cout<<"Esportazione terminata con successo per Cell2Ds.txt"<< endl;
	  	
	  	
	  if (! GeneraFileCell3D (mesh, outfilename3D))
	 {
		 cerr<<"Problemi di esportazione per Cell3Ds.txt"<< endl;
		 
	  }
	  else
	  	cout<<"Esportazione terminata con successo per Cell3Ds.txt"<< endl;
	  	
	  }
		 
/*-----------------------------------------------------------------------------------------------*/

/// Funzione che genera il file Cell0D

bool GeneraFileCell0D(const PolygonalMesh& mesh, const string& outfilename) {
    ofstream file(outfilename);
    if (!file.is_open()) return false;

    file << "Id, x, y, z, ShortPath\n";
    
    auto it = mesh.MarkerCell0Ds.find(1);
    list<unsigned int> ShortPath;
    if (it != mesh.MarkerCell0Ds.end()) {
        list<unsigned int> ShortPath = it->second;
        cout << "La chiave 1 è presente nella mappa." << endl;
        // Usa ShortPath come necessario
    } else {
        // La chiave 1 non è presente
        cout << "La chiave 1 non è presente nella mappa." << endl;
    }


    for (unsigned int i = 0; i < mesh.NumCell0Ds; ++i) {
        file << mesh.Cell0DsId[i] << ", "
             << setprecision(10) << mesh.Cell0DsCoordinates(0, i) << ", "
             << setprecision(10) << mesh.Cell0DsCoordinates(1, i) << ", "
			 << setprecision(10) << mesh.Cell0DsCoordinates(2, i) << ", ";
            if (find(ShortPath.begin(), ShortPath.end(), mesh.Cell0DsId[i]) != ShortPath.end()) {
                file << 1 << "\n";
            } else {
                file << 0 << "\n";
            }
    }
    file.close();
	
	
    return true;
}

/*-----------------------------------------------------------------------------------------------*/

/// Funzione che genera il file Cell1D

bool GeneraFileCell1D(const PolygonalMesh& mesh, const string& outfilename) {
    ofstream file(outfilename);
    if (!file.is_open()) {
        cerr << "Errore apertura file " << outfilename << endl;
        return false;
    }

    file << "Id, Id_Origine, Id_Fine, ShortPath\n";

    auto it = mesh.MarkerCell1Ds.find(1);
    list<unsigned int> ShortPath;
    if (it != mesh.MarkerCell1Ds.end()) {
        list<unsigned int> ShortPath = it->second;
        cout << "La chiave 1 è presente nella mappa." << endl;
        // Usa ShortPath come necessario
    } else {
        // La chiave 1 non è presente
        cout << "La chiave 1 non è presente nella mappa." << endl;
    }

    for (unsigned int i= 0; i < mesh.NumCell1Ds; ++i) {
        file << mesh.Cell1DsId[i] << ", " 
		<< setprecision(10) << mesh.Cell1DsExtrema(0,i)<< ", " 
		<<setprecision(10) << mesh.Cell1DsExtrema(1,i) << ", ";
        
        if (find(ShortPath.begin(), ShortPath.end(), mesh.Cell1DsId[i]) != ShortPath.end()) {
            file << 1 << "\n";
        } else {
            file << 0 << "\n";
        }
        
    }
    file.close();

    return true;
	
}

/*-----------------------------------------------------------------------------------------------*/

/// Funzione che genera il file Cell2D

bool GeneraFileCell2D(const PolygonalMesh& mesh, const string& outfilename) {
    ofstream file(outfilename);
    if (!file.is_open()) {
        cerr << "Errore apertura file " << outfilename << endl;
        return false;
    }

    file << "Id, NumVertici, NumSpigoli, Id_Vertici, Id_Spigoli\n";

    for (unsigned int i= 0; i < mesh.NumCell2Ds; ++i) {
        file << mesh.Cell2DsId[i] << ", "
		<< setprecision(10) << mesh.Cell2DsVertices[i].size() << ", "
        << setprecision(10) << mesh.Cell2DsEdges[i].size() << ", ";
		
		file << "[" ;
        // Vertici
        for (size_t j=0; j < mesh.Cell2DsVertices[i].size();j++)
            file << mesh.Cell2DsVertices[i][j] << "  ";
        // Spigoli
        file << "], [" ;
        for (size_t j=0; j<mesh.Cell2DsEdges[i].size();j++)
            file << mesh.Cell2DsEdges[i][j] << "  ";

        file << "]\n";
    }
    file.close();
    return true;
}

/*-----------------------------------------------------------------------------------------------*/

/// Funzione che genera il file Cell3D


bool GeneraFileCell3D(const PolygonalMesh& mesh, const string& outfilename) {
    ofstream file(outfilename);
    if (!file.is_open()) {
        cerr << "Errore apertura file " << outfilename << endl;
        return false;
    }

    file << "Id, NumVertici, NumSpigoli, NumFacce, Id_Vertici, Id_Spigoli, Id_Facce\n";
	
	for (unsigned int i = 0; i < mesh.NumCell3Ds; ++i) {
        file << mesh.Cell0DsId[i] << ", "
		<< setprecision(10) << mesh.Cell3DsVertices[i].size() << ", "
		<< setprecision(10) << mesh.Cell3DsEdges[i].size() << ", "
		<< setprecision(10) << mesh.Cell3DsFaces[i].size() << ", ";
		
		file << "[" ;
        // Vertici
        for (size_t j=0; j<mesh.Cell3DsVertices[i].size();j++)
            file << mesh.Cell3DsVertices[i][j] << "  ";
		
		// Spigoli
        file << "], [" ;
        for (size_t j=0; j<mesh.Cell3DsEdges[i].size();j++)
            file << mesh.Cell3DsEdges[i][j] << "  ";
		
		// Facce
        file << "], [" ;
        for (size_t j=0; j<mesh.Cell3DsFaces[i].size();j++)
            file << mesh.Cell3DsFaces[i][j] << "  ";

        file << "]\n";
		
    }
    file.close();
    return true;
}


/*---------------------------------/// DUALE /// -------------------------------------------------------------*/

/*
Function: calcola il CentroideCoordinate
input2: Cell0DsCoordinates
*/

static Vector3d CentroideCoordinate(const vector<int>& VerticesId, const MatrixXd& Cell0DsCoordinates) {
    Vector3d c(0.0, 0.0, 0.0);
    for (int vid : VerticesId) {   //itera su ogni punto (con coordinate in colonna)
        c += Cell0DsCoordinates.col(vid);
    }
    c /= static_cast<double>(VerticesId.size());   //media (statica_cast trasforma l'intero in double per evitare la diviosione intera)
    
	double norm = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);   // calcolo la norma per normalizzare
    Vector3d p(c[0]/norm, c[1]/norm, c[2]/norm);
	
    return p;     
    // VALUTARE round()
}

/*-----------------------------------------------------------------------------------------------*/
/*
Function: costruisce il duale
*/

void build_duale(const PolygonalMesh& geodetico, PolygonalMesh& Goldberg) {

    Goldberg = PolygonalMesh();   // serve per evitare una copia e lavorare direttamente per reference

    // === Cell0D: un vertice per ogni faccia geodetico ===
    map<unsigned int, unsigned int> mapGoldbergVerticiID;   //mappa per andare dalla faccia del geodetico al vertice del duale
    Goldberg.NumCell0Ds = geodetico.NumCell2Ds;   // proprietà del duale (scambio facce e vertici)
    Goldberg.Cell0DsId.resize(Goldberg.NumCell0Ds);   // adattare le dimensioni in base al NumCell0Ds
    Goldberg.Cell0DsCoordinates.resize(3, Goldberg.NumCell0Ds);

    for (unsigned int i = 0; i < geodetico.NumCell2Ds; ++i) {
        Goldberg.Cell0DsId[i] = i;   //Id nella mesh
        Vector3d baricentro = CentroideCoordinate(geodetico.Cell2DsVertices[i], geodetico.Cell0DsCoordinates);
        Goldberg.Cell0DsCoordinates.col(i) = baricentro;   //coordinate nella mesh
        mapGoldbergVerticiID[i] = i;   //mappatura
    }

    // === Cell2D: una faccia per ogni vertice geodetico ===
    vector<vector<int>> GoldbergFacce;   //per andare dal vertice del geodetico alla faccia di Goldberg

    // Mappa dei lati per trovare le facce che li condividono
    map<pair<int, int>, vector<int>> mapLatiFacce; // (lato) -> (facce che condividono questo lato)

    // Costruzione della mappa dei lati e delle facce adiacenti
    for (unsigned int fid = 0; fid < geodetico.NumCell2Ds; ++fid) {
        const auto& vertici = geodetico.Cell2DsVertices[fid];
        
        for (size_t i = 0; i < vertici.size(); ++i) {
            int a = vertici[i];
            int b = vertici[(i + 1) % vertici.size()];
            
            if (a > b) swap(a, b);  // Garantiamo che l'ordine sia (minore, maggiore)
            
            mapLatiFacce[make_pair(a, b)].push_back(fid); // Aggiungi faccia alla mappa
        }
    }

    // Ora costruisci le facce nel duale
    for (const auto& [lato, facce] : mapLatiFacce) {
        if (facce.size() == 2) {
            // Le due facce condividono questo lato, quindi collega i loro baricentri nel duale
            int fid1 = facce[0];
            int fid2 = facce[1];
            
            // Aggiungi i baricentri delle facce adiacenti come vertici nel duale
            int v1 = mapGoldbergVerticiID[fid1];
            int v2 = mapGoldbergVerticiID[fid2];
            
            // Crea uno spigolo tra i baricentri (collega le facce nel duale)
            GoldbergFacce.push_back({v1, v2}); // Qui creiamo una faccia composta da due vertici
        }
    }
    Goldberg.NumCell2Ds = static_cast<unsigned int>(GoldbergFacce.size());
    Goldberg.Cell2DsId.resize(Goldberg.NumCell2Ds);
    Goldberg.Cell2DsVertices = GoldbergFacce;
    Goldberg.Cell2DsNumEdges.resize(Goldberg.NumCell2Ds);
    Goldberg.Cell2DsEdges.resize(Goldberg.NumCell2Ds);

    // === Cell1D: crea spigoli univoci dalle facce ===
    map<pair<int, int>, int> mapEdgeGoldberg;
    int eid = 0;
    for (unsigned int fid = 0; fid < Goldberg.NumCell2Ds; ++fid) 
    {
        Goldberg.Cell2DsId[fid] = fid;   //mettiamo gli Id ad ogni faccia di Goldberg
        const auto& vertici = Goldberg.Cell2DsVertices[fid];
        Goldberg.Cell2DsNumEdges[fid] = vertici.size();   //numero di lati e vertici è lo stesso

        for (size_t i = 0; i < vertici.size(); ++i) {
            int a = vertici[i];
            int b = vertici[(i+1)%vertici.size()];
            if (a > b) swap(a, b);
            auto key = make_pair(a, b);
            if (!mapEdgeGoldberg.count(key)) {
                mapEdgeGoldberg[key] = eid;
                Goldberg.Cell1DsId.push_back(eid);
                ++eid;
            }
            Goldberg.Cell2DsEdges[fid].push_back(mapEdgeGoldberg[key]);   //mette ogni lato alla faccia corrispondente
        }
    }
    Goldberg.NumCell1Ds = Goldberg.Cell1DsId.size();
    Goldberg.Cell1DsExtrema.resize(2, Goldberg.NumCell1Ds);
    for (const auto& [key, id] : mapEdgeGoldberg) {
        Goldberg.Cell1DsExtrema(0, id) = key.first;
        Goldberg.Cell1DsExtrema(1, id) = key.second;
    }

    // === Cell3D: unica cella 3D contenente tutto ===
    Goldberg.NumCell3Ds = 1;
    Goldberg.Cell3DsId = {0};

    vector<int> vertici3D(Goldberg.NumCell0Ds);
    iota(vertici3D.begin(), vertici3D.end(), 0);
    vector<int> lati3D(Goldberg.NumCell1Ds);
    iota(lati3D.begin(), lati3D.end(), 0);
    vector<int> facce3D(Goldberg.NumCell2Ds);
    iota(facce3D.begin(), facce3D.end(), 0);

    Goldberg.Cell3DsVertices = {vertici3D};
    Goldberg.Cell3DsEdges = {lati3D};
    Goldberg.Cell3DsFaces = {facce3D};
} 

/*-----------------------------------------------------------------------------------------------*/

// Funzione per calcolare il cammino minimo

void Dijkstra(PolygonalMesh& mesh, unsigned int vertice_iniziale, unsigned int vertice_finale)   //sono gli id dei vertici 
{
    // Mappa di adiacenza per il grafo: key = id del vertice, value = lista di id dei vertici adiacenti
    map<unsigned int, list<unsigned int>> ListaVerticiAdiacenti;
    
    // Creiamo la lista di adiacenza dai dati in Cell2DsVertices (cellette 2D)
    for (unsigned int i = 0; i < mesh.NumCell2Ds; ++i) {
        const auto& vertices = mesh.Cell2DsVertices[i];
        for (size_t j = 0; j < vertices.size(); ++j) {
            unsigned int v1 = vertices[j];
            unsigned int v2 = vertices[(j + 1) % vertices.size()]; // Arco da v1 a v2 (ciclo per il bordo)
            
            ListaVerticiAdiacenti[v1].push_back(v2);
            ListaVerticiAdiacenti[v2].push_back(v1); // Grafo non orientato
        }
    }

    // Algoritmo di Dijkstra per trovare il cammino minimo

    map<unsigned int, double> dist;   //vertici e relative distanze rispetto al nodo sorgente (vertice_iniziale)
    map<unsigned int, unsigned int> pred;   //vertici e relativi nodi precedenti lungo il cammino minimo
    priority_queue< pair<double, unsigned int>, vector<pair<double, unsigned int>>, greater<> > pq;
    
    // Inizializza la distanza del vertice di partenza
    for (auto& entrata : ListaVerticiAdiacenti) {
        dist[entrata.first] = numeric_limits<double>::infinity();   //serve per impostare un infinito
    }
    dist[vertice_iniziale] = 0;
    pq.push({0, vertice_iniziale});

    for (size_t i=0; i<dist.size(); i++){
        pq.push({dist[i], i});
    }
    
    while (!pq.empty()) {
        unsigned int u = pq.top().second;   //u=vertice
        double d = pq.top().first;   //d=distanza di u dal vertice_iniziale
        pq.pop();   //tolgo l'elemento dalla coda con priorità
        
        // Se la distanza per il vertice corrente è maggiore della distanza già trovata, salta
        if (d > dist[u]) continue;   //se soddisfatta, continuo col for
        
        // Esamina i vicini
        for (unsigned int v : ListaVerticiAdiacenti[u]) {
            double nuova_dist = dist[u] + 1; // Peso dell'arco (1 per semplicità, puoi usare distanze Euclidee se le hai)
            if (nuova_dist < dist[v]) {
                dist[v] = nuova_dist;
                pred[v] = u;
                pq.push({nuova_dist, v});
            }
        }
    }

    // Ricostruisci il cammino minimo partendo da vertice2
    vector<unsigned int> path;
    for (unsigned int u = vertice_finale; u != vertice_iniziale; u = pred[u]) {
        if (pred.find(u) == pred.end()) {
            cout << "Nessun cammino trovato!" << endl;
            return;
        }
        path.push_back(u);
    }
    path.push_back(vertice_iniziale);
    reverse(path.begin(), path.end());


    // Marca i vertici nel cammino
    /*for (unsigned int v : path) {
        // Verifica prima se esiste la chiave 1 nella mappa
        cout << "Verifica MarkerCell1Ds prima dell'inserimento: ";
        for (auto& entry : mesh.MarkerCell1Ds) {
            cout << entry.first << ": ";
            for (unsigned int marker : entry.second) {
                cout << marker << " ";
            }
            cout << endl;
        }
    }*/
    
    for (unsigned int v : path) {
            mesh.MarkerCell0Ds[1].push_back(v);
            cout << "vertice" << v <<endl;
        }

    /*cout << "Verifica MarkerCell1Ds dopo l'inserimento: ";
        for (auto& entry : mesh.MarkerCell1Ds) {
            cout << entry.first << ": ";
            for (unsigned int marker : entry.second) {
                cout << marker << " ";
            }
            cout << endl;*/

        /*const auto it = mesh.MarkerCell1Ds.find(1);   // uso "auto" poichè it è un iteratore che non so dove punta 
            if(it == mesh.MarkerCell1Ds.end())   // aggiungo un nuovo marker
            {
                mesh.MarkerCell1Ds.insert({1, {v}});
                cout << "vertice" << v <<endl;
            }
            else   // aggiungo un elemento al marker
            {
                mesh.MarkerCell1Ds[1].push_back(v);
                cout << "vertice" << v <<endl;
            }
        }*/

    // Marca gli archi nel cammino
    for (size_t i = 1; i < path.size(); ++i) {
        unsigned int u = path[i-1];
        unsigned int v = path[i];
        // Trova l'arco tra u e v (Cell1D)
        for (unsigned int j = 0; j < mesh.NumCell1Ds; ++j) {
            if ((mesh.Cell1DsExtrema(0, j) == u && mesh.Cell1DsExtrema(1, j) == v) ||
                (mesh.Cell1DsExtrema(0, j) == v && mesh.Cell1DsExtrema(1, j) == u)) {
                const auto it = mesh.MarkerCell1Ds.find(1);   // uso "auto" poichè it è un iteratore che non so dove punta 
                if(it == mesh.MarkerCell1Ds.end())   // aggiungo un nuovo marker
                {
                    mesh.MarkerCell1Ds.insert({1, {mesh.Cell1DsId[j]}});
                    cout << "spigolo" << j <<endl;
                }
                else   // aggiungo un elemento al marker
                {
                    mesh.MarkerCell1Ds[1].push_back(mesh.Cell1DsId[j]);
                    cout << "spigolo" << j <<endl;
                }
                break;
            }
        }
    }
}





