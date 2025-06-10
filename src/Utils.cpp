#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
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
#include <set>
#include <algorithm>

using namespace PolyhedralLibrary;
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
                           PolyhedralMesh& mesh) 
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

void build_classe_1(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolyhedralMesh& mesh) {
    
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

/*
// Funzione : inserisci_vertice : normalizza il vertice , assegna l'id , inserisce nella mesh l'id e le coordinate nel mesh
// Inpputs :
coord : coordinate del vertice
mesh
// Outputs : id del vertice
*/
unsigned int inserisci_vertice(array<double,3> coord, 
                               PolyhedralMesh& mesh) {
    array<double,3> norm = Normalizza(coord);
    for (unsigned int i = 0; i < mesh.NumCell0Ds; ++i) {
        if ((abs(mesh.Cell0DsCoordinates(0,i) - norm[0]) < 1e-8) &&   //controllo se è già nel Cell0DsCoordinates
            (abs(mesh.Cell0DsCoordinates(1,i) - norm[1]) < 1e-8) &&
            (abs(mesh.Cell0DsCoordinates(2,i) - norm[2]) < 1e-8)) { return i;}
    
    unsigned int id = mesh.NumCell0Ds++;
    mesh.Cell0DsId.push_back(id);   // Cell0DsID
    mesh.Cell0DsCoordinates.conservativeResize(3, id+1);   // Cell0DsCoordinates
    mesh.Cell0DsCoordinates(0,id) = norm[0];
    mesh.Cell0DsCoordinates(1,id) = norm[1];
    mesh.Cell0DsCoordinates(2,id) = norm[2];
    return id;
    }
}

/*-----------------------------------------------------------------------------------------------*/

/*
// Funzione : Baricentro : calcola in coordinate il baricentro di un triangolo
// Inputs : A,B,C : vertici del triangolo (in coord)
// Outputs : coordinate del baricentro di quel triangolo
*/
array<double,3> Baricentro(const array<double,3>& A, const array<double,3>& B, const array<double,3>& C) {
return {(A[0]+B[0]+C[0])/3.0, (A[1]+B[1]+C[1])/3.0, (A[2]+B[2]+C[2])/3.0};
}

/*-----------------------------------------------------------------------------------------------*/
/*
//Funzione : PuntoMedio : calcola in coordinate  il punto medio di un lato
// Inputs : A,B : estremi di un lato (in coord)
// Outputs : punto medio in coordinate
*/
array<double,3> PuntoMedio(const array<double,3>& A, const array<double,3>& B) {
return {(A[0]+B[0])/2.0, (A[1]+B[1])/2.0, (A[2]+B[2])/2.0};
}

/*-----------------------------------------------------------------------------------------------*/

/* 
// Funzione : mesh secondo la triangolazione di classe 2 con b = c con b>0 (classe II), considerando:
// Inputs : p, q, b, c==b con b > 0
// Outputs : riempimento del mesh in PolyhedralMesh
*/
void build_classe_2(unsigned int p, 
                    unsigned int q, 
                    unsigned int b, 
                    unsigned int c, 
                    PolyhedralMesh& mesh) {

    map<array<int,3>, unsigned int> coordBari_id;   // key: vettori per coord bari -> id del bari  
    map<pair<unsigned int,unsigned int>, unsigned int> archi_mappa;   // key: coppia di id di vertici -> id dell'arco
    set<pair<unsigned int,unsigned int>> archi;   // raccolta di tutti gli archi identificati dagli id dei vertici

    vector<array<double, 3>> vertici_iniziali;
    vector<array<int, 3>> facce_iniziali;

    if (q == 5) build_ico(vertici_iniziali, facce_iniziali);
    else if (q == 4) build_octa(vertici_iniziali, facce_iniziali);
    else build_tetra(vertici_iniziali, facce_iniziali);

    for (const auto& face : facce_iniziali) {
        const auto& A = vertici_iniziali[face[0]];
        const auto& B = vertici_iniziali[face[1]];
        const auto& C = vertici_iniziali[face[2]]; 

        /// "GRIGLIA" COME NEL build_classe_1
        for (int i = 0; i <= (int)b; ++i) {
            for (int j = 0; j <= (int)(b - i); ++j) {
                int k = b - i - j;   // regola base delle coord baricentrice
                double x = (i * vertici_iniziali[0][0] + j * vertici_iniziali[1][0] + k * vertici_iniziali[2][0]) / double(b);   //combinaz convessa del triangolo iniziale
                double y = (i * vertici_iniziali[0][1] + j * vertici_iniziali[1][1] + k * vertici_iniziali[2][1]) / double(b);
                double z = (i * vertici_iniziali[0][2] + j * vertici_iniziali[1][2] + k * vertici_iniziali[2][2]) / double(b);
                array<double,3> P = {x,y,z};
                coordBari_id[{i,j,k}] = inserisci_vertice(P, mesh);   // assegno un id alle coord
            }
        }

        vector<unsigned int> bari_id;   // insieme di id dei baricentri

        vector<tuple<unsigned int,unsigned int,unsigned int>> sotto_triangoli;   // vettore di sottotriangoli rappresentati dagli id
        map<pair<unsigned int,unsigned int>, unsigned int> punti_medi;    // key: lato con estremi vertici -> id

        /// SOTTO-TRIANGOLI
        for (int i = 0; i < (int)b; ++i) {
            for (int j = 0; j < (int)(b - i); ++j) {
                int k = b - i - j;
                unsigned int a = coordBari_id[{i,j,k}];
                unsigned int b1 = coordBari_id[{i+1,j,k-1}];
                unsigned int c = coordBari_id[{i,j+1,k-1}];
                sotto_triangoli.emplace_back(a,b1,c);
                if (i + j < (int)(b - 1)) {
                    unsigned int d = coordBari_id[{i+1,j+1,k-2}];
                    sotto_triangoli.emplace_back(b1,d,c);
                }
            }
        }

        /// LATI CONDIVISI DA DUE SOTTOTRIANGOLI
        map<pair<unsigned int,unsigned int>, int> lati_condivisi;   // key: lato -> ripetizioni

        for (const auto& tri : sotto_triangoli) {
            array<unsigned int, 3> v = {get<0>(tri), get<1>(tri), get<2>(tri)};
            for (int i = 0; i < 3; ++i) {
                unsigned int a = v[i];
                unsigned int b = v[(i+1)%3];
                if (a > b) swap(a, b);
                lati_condivisi[{a, b}]++;
            }
        }

        /// ANALISI DI OGNI SOTTO-TRIANGOLO
        for (auto& tri : sotto_triangoli) {   
            unsigned int i0 = get<0>(tri);   // id del primo vertice
            unsigned int i1 = get<1>(tri);
            unsigned int i2 = get<2>(tri);

            /// CALCOLO DEL BARICENTRO
            array<double,3> A = {mesh.Cell0DsCoordinates(0,i0), mesh.Cell0DsCoordinates(1,i0), mesh.Cell0DsCoordinates(2,i0)};
            array<double,3> B = {mesh.Cell0DsCoordinates(0,i1), mesh.Cell0DsCoordinates(1,i1), mesh.Cell0DsCoordinates(2,i1)};
            array<double,3> C = {mesh.Cell0DsCoordinates(0,i2), mesh.Cell0DsCoordinates(1,i2), mesh.Cell0DsCoordinates(2,i2)};
            array<double,3> G = Baricentro(A,B,C);   // calcolo bari del sottotriangolo
            unsigned int G_id = inserisci_vertice(G, mesh);   // Cell0DsCoordinates e Cell0DsId del baricentro
            bari_id.push_back(G_id);

            /// LATI COI PUNTI_MEDI
            array<pair<unsigned int,unsigned int>,3> lati = {make_pair(i0,i1),make_pair(i1,i2),make_pair(i2,i0)};   //lati del sottotriangolo
            vector<int> lati_id;   // id lati
            for (auto [u,v] : lati) {
                if (u > v) swap(u,v);   //  ordine (min, max)

                /// PUNTI MEDI
                if (lati_condivisi[{u,v}] == 1) {
                    array<double,3> P1 = {mesh.Cell0DsCoordinates(0,u), mesh.Cell0DsCoordinates(1,u), mesh.Cell0DsCoordinates(2,u)};
                    array<double,3> P2 = {mesh.Cell0DsCoordinates(0,v), mesh.Cell0DsCoordinates(1,v), mesh.Cell0DsCoordinates(2,v)};
                    array<double,3> M = PuntoMedio(P1, P2);
                    unsigned int M_id = inserisci_vertice(M,mesh);

                    /// LATI BARI-PUNTO_MEDIO
                    pair<unsigned int,unsigned int> BM_id = {min(M_id,G_id), max(M_id,G_id)};   
                    if (!archi_mappa.count(BM_id)) {
                        unsigned int id = archi_mappa.size();
                        archi_mappa[BM_id] = id;
                        archi.insert(BM_id);
                    }

                    /// LATI PUNTO_MEDIO-VERTICI_GENERATORI
                    pair<unsigned int,unsigned int> uG_id = {min(u,G_id), max(u,G_id)}; 
                    if (!archi_mappa.count(uG_id)) {
                        unsigned int id = archi_mappa.size();
                        archi_mappa[uG_id] = id;
                        archi.insert(uG_id);
                    }
                    pair<unsigned int,unsigned int> vG_id = {min(v,G_id), max(v,G_id)};
                    if (!archi_mappa.count(vG_id)) {
                        unsigned int id = archi_mappa.size();
                        archi_mappa[vG_id] = id;
                        archi.insert(vG_id);
                    }
                }
            }

            // LATI BARI-VERTICI_SOTTOTRIANGOLO
            for (auto v : {i0, i1, i2}) {
                pair<unsigned int,unsigned int> BV_id = {min(G_id,v), max(G_id,v)};
                if (!archi_mappa.count(BV_id)) {
                    unsigned int id = archi_mappa.size();
                    archi_mappa[BV_id] = id;
                    archi.insert(BV_id);
                }
            }
        }

        /// LATI BARI-BARI ADIACENTI
        for (size_t i = 0; i < sotto_triangoli.size(); ++i) {
            for (size_t j = i+1; j < sotto_triangoli.size(); ++j) {
                set<unsigned int> v1 = {get<0>(sotto_triangoli[i]), get<1>(sotto_triangoli[i]), get<2>(sotto_triangoli[i])};
                set<unsigned int> v2 = {get<0>(sotto_triangoli[j]), get<1>(sotto_triangoli[j]), get<2>(sotto_triangoli[j])};
                vector<unsigned int> vertici_condivisi;
                set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(vertici_condivisi));   //i vertici identici
                if (vertici_condivisi.size() == 2) {   // allora c'è  un lato in comune
                    unsigned int a = bari_id[i];
                    unsigned int b = bari_id[j];
                    pair<unsigned int,unsigned int> key = {min(a,b), max(a,b)};
                    if (!archi_mappa.count(key)) {
                        unsigned int id = archi_mappa.size();
                        archi_mappa[key] = id;
                        archi.insert(key);
                    }
                }
            }
        }

        /// LATI NELLA MESH
        mesh.Cell1DsExtrema.resize(2, archi.size());
        map<pair<int,int>, int> mappa_id_archi;
        mesh.NumCell1Ds = 0;   // NumCell1Ds
        for (const auto& [a,b] : archi) {
            mesh.Cell1DsId.push_back(mesh.NumCell1Ds);   // Cell1DsId
            mesh.Cell1DsExtrema(0, mesh.NumCell1Ds) = a;   // Cell1DsExtrema
            mesh.Cell1DsExtrema(1, mesh.NumCell1Ds) = b;
            mappa_id_archi[make_pair(a,b)] = mesh.NumCell1Ds;
            mesh.NumCell1Ds++;
        }
        
        /// FACCE NELLA MESH
        mesh.NumCell2Ds = 0;
        for(const auto& [a,b] : archi) {
            unsigned int u = a;
            unsigned int v = b;

            vector<pair<unsigned int, unsigned int>> lati_consecutivi;   // raccolta di lati che hanno come un estremo v
            for(const auto& [c,d] : archi) {
                if (c==v || d==v && c!=u) {
                    lati_consecutivi.push_back(make_pair(c,d));

                    for (const auto& [i,j] : lati_consecutivi) {
                        if (i==u || j==u) { 
                            vector<int> triangolo = {u,v,i,j};   // faccia
                            sort(triangolo.begin(), triangolo.end());
                            auto it = unique(triangolo.begin(), triangolo.end());   //individuo i vertici ripetuti
                            triangolo.erase(it, triangolo.end());   // elimino i vertici ripetuti

                            pair<int,int> lato1  = make_pair(triangolo[0], triangolo[1]) ;
                            pair<int,int> lato2  = make_pair(triangolo[1], triangolo[2]) ;
                            pair<int,int> lato3  = make_pair(triangolo[0], triangolo[2]) ;

                            if (find(mesh.Cell2DsVertices.begin(), mesh.Cell2DsVertices.end(), triangolo)== mesh.Cell2DsVertices.end() ){
                                mesh.Cell2DsId.push_back(mesh.NumCell2Ds);   // Cell2DsId
                                mesh.Cell2DsVertices.push_back(triangolo);   // Cell2DsVertices
                                mesh.Cell2DsEdges.push_back({mappa_id_archi[lato1],mappa_id_archi[lato2], mappa_id_archi[lato3]});   //Cell2DsEdges
                                mesh.Cell2DsNumEdges.push_back(3);   // Cell2DsNumEdges
                                mesh.NumCell2Ds++;
                            }
                        }
                    }
                }
            }
        }
    }

	// CELL3D (UNICA)
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
void GeneraTuttiFile(const PolyhedralMesh& mesh, 
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

bool GeneraFileCell0D(const PolyhedralMesh& mesh, const string& outfilename) {
    ofstream file(outfilename);
    if (!file.is_open()) return false;

    file << "Id, x, y, z\n";
    
    // CONTROLLO MARKER (da eliminare)
    auto it = mesh.MarkerCell0Ds.find(1);
    if (it != mesh.MarkerCell0Ds.end()) {
        cout << "La chiave 1 è presente nella mappa." << endl;
    } else {
        cout << "La chiave 1 non è presente nella mappa." << endl;
    }


    for (unsigned int i = 0; i < mesh.NumCell0Ds; ++i) {   //stampa sul outputFile
        file << mesh.Cell0DsId[i] << ", "
             << setprecision(10) << mesh.Cell0DsCoordinates(0, i) << ", "
             << setprecision(10) << mesh.Cell0DsCoordinates(1, i) << ", "
			 << setprecision(10) << mesh.Cell0DsCoordinates(2, i) << "\n";
            
    }

    file.close();
	
	
    return true;
}

/*-----------------------------------------------------------------------------------------------*/

/// Funzione che genera il file Cell1D

bool GeneraFileCell1D(const PolyhedralMesh& mesh, const string& outfilename) {
    ofstream file(outfilename);
    if (!file.is_open()) {
        cerr << "Errore apertura file " << outfilename << endl;
        return false;
    }

    file << "Id, Id_Origine, Id_Fine\n";

    // CONTROLLO MARKER (da eliminare)
    auto it = mesh.MarkerCell1Ds.find(1);
    if (it != mesh.MarkerCell1Ds.end()) {
        cout << "La chiave 1 è presente nella mappa." << endl;
    } else {
        cout << "La chiave 1 non è presente nella mappa." << endl;
    }

    for (unsigned int i= 0; i < mesh.NumCell1Ds; ++i) {   // stampa sul outputFile
        file << mesh.Cell1DsId[i] << ", " 
		<< setprecision(10) << mesh.Cell1DsExtrema(0,i)<< ", " 
		<<setprecision(10) << mesh.Cell1DsExtrema(1,i) << "\n";
    }
    file.close();

    return true;
	
}

/*-----------------------------------------------------------------------------------------------*/

/// Funzione che genera il file Cell2D

bool GeneraFileCell2D(const PolyhedralMesh& mesh, const string& outfilename) {
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


bool GeneraFileCell3D(const PolyhedralMesh& mesh, const string& outfilename) {
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

void build_duale(const PolyhedralMesh& geodetico, PolyhedralMesh& Goldberg) {

    Goldberg = PolyhedralMesh();   // serve per evitare una copia e lavorare direttamente per reference

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

void Dijkstra(PolyhedralMesh& mesh, unsigned int vertice_iniziale, unsigned int vertice_finale)   //sono gli id dei vertici 
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
            double dist_uv = sqrt(pow(mesh.Cell0DsCoordinates(0,u)-mesh.Cell0DsCoordinates(0,v),2) +   //distanza euclidea uv
                                  pow(mesh.Cell0DsCoordinates(1,u)-mesh.Cell0DsCoordinates(1,v),2) +
                                  pow(mesh.Cell0DsCoordinates(2,u)-mesh.Cell0DsCoordinates(2,v),2) );
            double nuova_dist = dist[u] + dist_uv;   // dist_u0 + dist_uv
            if (nuova_dist < dist[v]) {
                dist[v] = nuova_dist;
                pred[v] = u;
                pq.push({nuova_dist, v});
            }
        }
    }

    // Ricostruisci il cammino minimo partendo da vertice_finale
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

    //memorizzo il numero degli archi del path e la somma totale delle distanza
    mesh.num_archiPath = path.size()-1;   //poichè path contiene vertici
    mesh.lunghezza_Path = dist[vertice_finale];
    
    // marcare vertici del cammino minimo
    for (unsigned int v : path) {
            mesh.MarkerCell0Ds[1].push_back(v);
            cout << "vertice" << v <<endl;
        }

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
