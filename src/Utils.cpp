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
#include <unordered_set>

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
                           PolyhedralMesh& mesh) {
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

pair<int, bool> get_or_add_edge(
    unsigned int a, unsigned int b,
    map<pair<unsigned int, unsigned int>, unsigned int>& map1D,
    PolyhedralMesh& mesh) {
    pair<unsigned int, unsigned int> key = (a < b) ? make_pair(a, b) : make_pair(b, a);
    bool orientation = (a < b); // true se orientamento coerente

    if (map1D.count(key)) {
        int eid = map1D[key];
        bool actual_orientation = (a == mesh.Cell1DsExtrema(0, eid));
        return {eid, actual_orientation};
    }

    // Aggiunta nuovo lato
    int eid = mesh.Cell1DsId.size();
    map1D[key] = eid;
    mesh.Cell1DsId.push_back(eid);
    mesh.Cell1DsExtrema.conservativeResize(2, eid + 1);
    mesh.Cell1DsExtrema(0, eid) = key.first;
    mesh.Cell1DsExtrema(1, eid) = key.second;

    return {eid, orientation};
}

/*-----------------------------------------------------------------------------------------------*/

void build_classe_1(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolyhedralMesh& mesh) {
    using namespace std;

    vector<array<double, 3>> vertici_iniziali;
    vector<array<int, 3>> facce_iniziali;

    if (q == 5) build_ico(vertici_iniziali, facce_iniziali);
    else if (q == 4) build_octa(vertici_iniziali, facce_iniziali);
    else build_tetra(vertici_iniziali, facce_iniziali);

    map<array<double, 3>, unsigned int> map0D;
    map<pair<unsigned int, unsigned int>, unsigned int> map1D;
    unsigned int fid = 0;

    for (const auto& face : facce_iniziali) {
        const auto& A = vertici_iniziali[face[0]];
        const auto& B = vertici_iniziali[face[1]];
        const auto& C = vertici_iniziali[face[2]];

        vector<vector<unsigned int>> griglia(b + 1);
        for (unsigned int i = 0; i <= b; ++i) {
            griglia[i].resize(i + 1);
            for (unsigned int j = 0; j <= i; ++j) {
                double u = 1.0 - static_cast<double>(i) / b;
                double v = static_cast<double>(i - j) / b;
                double w = static_cast<double>(j) / b;
                array<double, 3> P = {
                    u * A[0] + v * B[0] + w * C[0],
                    u * A[1] + v * B[1] + w * C[1],
                    u * A[2] + v * B[2] + w * C[2]
                };
                griglia[i][j] = salva_vertice_norm(P, map0D, mesh);
            }
        }

        for (unsigned int i = 0; i < b; ++i) {
            for (unsigned int j = 0; j < i + 1; ++j) {
                unsigned int v1 = griglia[i][j];
                unsigned int v2 = griglia[i + 1][j];
                unsigned int v3 = griglia[i + 1][j + 1];

                mesh.Cell2DsId.push_back(fid++);
                mesh.Cell2DsVertices.push_back({(int)v1, (int)v2, (int)v3});

                vector<int> edge_ids;

                for (int k = 0; k < 3; ++k) {
                    unsigned int a = mesh.Cell2DsVertices.back()[k];
                    unsigned int b = mesh.Cell2DsVertices.back()[(k + 1) % 3];
                    auto lato = make_pair(a, b);
                    auto lato_inverso = make_pair(b, a);
                    
                    if (map1D.count(lato) == 0) {
                        if (map1D.count(lato_inverso) == 0) {   // se il lato, nè il suo inverso, è presente -> ce lo aggiungo
                            map1D[lato] = mesh.Cell1DsId.size();
                            mesh.Cell1DsId.push_back(map1D[lato]);
                            mesh.Cell1DsExtrema.conservativeResize(2, map1D[lato] + 1);
                            mesh.Cell1DsExtrema(0, map1D[lato]) = lato.first;
                            mesh.Cell1DsExtrema(1, map1D[lato]) = lato.second; 
                            edge_ids.push_back(map1D[lato]);
                        } else if (map1D.count(lato_inverso) != 0) {
                            edge_ids.push_back(map1D[lato_inverso]);
                        }
                    } else if(map1D.count(lato) != 0) {
                        edge_ids.push_back(map1D[lato]);
                    }
                }

                mesh.Cell2DsEdges.push_back(edge_ids);
                mesh.Cell2DsNumEdges.push_back(3);

                if (j < i) {
                    unsigned int v4 = griglia[i][j + 1];
                    mesh.Cell2DsId.push_back(fid++);
                    mesh.Cell2DsVertices.push_back({(int)v1, (int)v3, (int)v4});

                    vector<int> edge_ids2;

                    for (int k = 0; k < 3; ++k) {
                    unsigned int a = mesh.Cell2DsVertices.back()[k];
                    unsigned int b = mesh.Cell2DsVertices.back()[(k + 1) % 3];
                    auto lato = make_pair(a, b);
                    auto lato_inverso = make_pair(b, a);
                    
                    if (map1D.count(lato) == 0) {
                        if (map1D.count(lato_inverso) == 0) {   // se il lato, nè il suo inverso, è presente -> ce lo aggiungo
                            map1D[lato] = mesh.Cell1DsId.size();
                            mesh.Cell1DsId.push_back(map1D[lato]);
                            mesh.Cell1DsExtrema.conservativeResize(2, map1D[lato] + 1);
                            mesh.Cell1DsExtrema(0, map1D[lato]) = lato.first;
                            mesh.Cell1DsExtrema(1, map1D[lato]) = lato.second;
                            edge_ids2.push_back(map1D[lato]);
                        } else if (map1D.count(lato_inverso) != 0) {
                            edge_ids2.push_back(map1D[lato_inverso]);
                        }
                    } else if(map1D.count(lato) != 0) {
                        edge_ids2.push_back(map1D[lato]);
                    }
                }

                mesh.Cell2DsEdges.push_back(edge_ids2);
                mesh.Cell2DsNumEdges.push_back(3);
                }
            }
        }
    }

    // Post-processing: Cell3D
    mesh.NumCell0Ds = mesh.Cell0DsId.size();
    mesh.NumCell1Ds = mesh.Cell1DsId.size();
    mesh.NumCell2Ds = mesh.Cell2DsId.size();

    mesh.NumCell3Ds = 1;
    mesh.Cell3DsId.push_back(0);

    vector<int> cell3D_vertices(mesh.NumCell0Ds);
    iota(cell3D_vertices.begin(), cell3D_vertices.end(), 0);
    mesh.Cell3DsVertices.push_back(cell3D_vertices);

    vector<int> cell3D_edges(mesh.NumCell1Ds);
    iota(cell3D_edges.begin(), cell3D_edges.end(), 0);
    mesh.Cell3DsEdges.push_back(cell3D_edges);

    vector<int> cell3D_faces(mesh.NumCell2Ds);
    iota(cell3D_faces.begin(), cell3D_faces.end(), 0);
    mesh.Cell3DsFaces.push_back(cell3D_faces);

    // Ottimizzazione memoria
    mesh.Cell2DsId.shrink_to_fit();
    mesh.Cell2DsEdges.shrink_to_fit();
    mesh.Cell2DsVertices.shrink_to_fit();
    mesh.Cell2DsNumEdges.shrink_to_fit();
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
    
    vector<array<double, 3>> vertici_iniziali;
    vector<array<int, 3>> facce_iniziali;

    if (q == 5) build_ico(vertici_iniziali, facce_iniziali);
    else if (q == 4) build_octa(vertici_iniziali, facce_iniziali);
    else build_tetra(vertici_iniziali, facce_iniziali);

    map<pair<unsigned int,unsigned int>, unsigned int> map1D_2;   // key: coppia di id di vertici -> id dell'arco
    map<array<double, 3>, unsigned int> map0D_2;
    unsigned int fid = 0;

    
    set<pair<unsigned int,unsigned int>> archi;   // raccolta di tutti gli archi identificati dagli id dei vertici

    for (const auto& face : facce_iniziali) {
        const auto& A = vertici_iniziali[face[0]];
        const auto& B = vertici_iniziali[face[1]];
        const auto& C = vertici_iniziali[face[2]]; 

        map<array<int,3>, unsigned int> coordBari_id;   // key: vettori per coord bari -> id del bari  
        coordBari_id.clear();

        /// "GRIGLIA" COME NEL build_classe_1
        for (int i = 0; i <= (int)b; ++i) {
            for (int j = 0; j <= (int)(b - i); ++j) {
                int k = b - i - j;   // regola base delle coord baricentrice
                double x = (i * A[0] + j * B[0] + k * C[0]) / double(b);   //combinaz convessa del triangolo iniziale
                double y = (i * A[1] + j * B[1] + k * C[1]) / double(b);
                double z = (i * A[2] + j * B[2] + k * C[2]) / double(b);
                array<double,3> P = {x,y,z};
                coordBari_id[{i,j,k}] = salva_vertice_norm(P, map0D_2, mesh); // assegno un id alle coord
            }
        }


        vector<tuple<unsigned int,unsigned int,unsigned int>> sotto_triangoli;   // vettore di sottotriangoli rappresentati dagli id

        /// SOTTO-TRIANGOLI
        for (int i = 0; i < (int)b; ++i) {
            for (int j = 0; j < (int)(b - i); ++j) {
                int k = b - i - j;
                unsigned int v1 = coordBari_id[{i,j,k}];
                unsigned int v2 = coordBari_id[{i+1,j,k-1}];
                unsigned int v3 = coordBari_id[{i,j+1,k-1}];
                sotto_triangoli.emplace_back(v1, v2, v3);
                if (i + j < (int)(b - 1)) {
                    unsigned int v4 = coordBari_id[{i+1,j+1,k-2}];
                    sotto_triangoli.emplace_back(v2, v4, v3);
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
                lati_condivisi[{a, b}]++;   //incremento le volte in cui incontro l'arco a-b
            }
        }

        vector<unsigned int> bari_id;   // insieme di id dei baricentri

        /// ANALISI DI OGNI SOTTO-TRIANGOLO
        for (auto& tri : sotto_triangoli) {   
            unsigned int v1 = get<0>(tri);   // id del primo vertice
            unsigned int v2 = get<1>(tri);
            unsigned int v3 = get<2>(tri);

            /// CALCOLO DEL BARICENTRO
            array<double,3> A = {mesh.Cell0DsCoordinates(0,v1), mesh.Cell0DsCoordinates(1,v1), mesh.Cell0DsCoordinates(2,v1)};
            array<double,3> B = {mesh.Cell0DsCoordinates(0,v2), mesh.Cell0DsCoordinates(1,v2), mesh.Cell0DsCoordinates(2,v2)};
            array<double,3> C = {mesh.Cell0DsCoordinates(0,v3), mesh.Cell0DsCoordinates(1,v3), mesh.Cell0DsCoordinates(2,v3)};
            array<double,3> G = Baricentro(A,B,C);   // calcolo bari del sottotriangolo
            unsigned int G_id = salva_vertice_norm(G, map0D_2, mesh);   // Cell0DsCoordinates e Cell0DsId del baricentro
            bari_id.push_back(G_id);

            // LATI BARI-VERTICI_SOTTOTRIANGOLO
            for (auto v : {v1, v2, v3}) {
                unsigned int a = min(G_id,v);
                unsigned int b = max(G_id,v);
                pair<unsigned int,unsigned int> GV_key = {a,b};
                if (!map1D_2.count(GV_key)) {
                    map1D_2[GV_key] = mesh.Cell1DsId.size();
                    mesh.Cell1DsId.push_back(map1D_2[GV_key]);
                    mesh.Cell1DsExtrema.conservativeResize(2, map1D_2[GV_key] + 1);
                    mesh.Cell1DsExtrema(0, map1D_2[GV_key]) = a;
                    mesh.Cell1DsExtrema(1, map1D_2[GV_key]) = b;

                    archi.insert(GV_key);
                }
            }

            /// LATI COI PUNTI_MEDI
            array<pair<unsigned int,unsigned int>,3> lati = {make_pair(v1,v2),make_pair(v2,v3),make_pair(v3,v1)};   //lati del sottotriangolo
        
            for (auto [u,v] : lati) {
                if (u > v) swap(u,v);   //  ordine (min, max)

                /// PUNTI MEDI
                if (lati_condivisi[{u,v}] == 1) {
                    array<double,3> P1 = {mesh.Cell0DsCoordinates(0,u), mesh.Cell0DsCoordinates(1,u), mesh.Cell0DsCoordinates(2,u)};
                    array<double,3> P2 = {mesh.Cell0DsCoordinates(0,v), mesh.Cell0DsCoordinates(1,v), mesh.Cell0DsCoordinates(2,v)};
                    array<double,3> M = PuntoMedio(P1, P2);
                    unsigned int M_id = salva_vertice_norm(M, map0D_2, mesh);

                    /// LATI BARI-PUNTO_MEDIO
                    unsigned int a = min(M_id,G_id);
                    unsigned int b = max(M_id,G_id);
                    pair<unsigned int,unsigned int> MG_key = {a,b};
                    if (!map1D_2.count(MG_key)) {
                        map1D_2[MG_key] = mesh.Cell1DsId.size();
                        mesh.Cell1DsId.push_back(map1D_2[MG_key]);
                        mesh.Cell1DsExtrema.conservativeResize(2, map1D_2[MG_key] + 1);
                        mesh.Cell1DsExtrema(0, map1D_2[MG_key]) = a;
                        mesh.Cell1DsExtrema(1, map1D_2[MG_key]) = b;

                        archi.insert(MG_key);
                    }

                    /// LATI PUNTO_MEDIO-VERTICI_GENERATORI
                    unsigned int c = min(u,M_id);
                    unsigned int d = max(u,M_id);
                    pair<unsigned int,unsigned int> uM_key = {c,d};
                    if (!map1D_2.count(uM_key)) {
                        map1D_2[uM_key] = mesh.Cell1DsId.size();
                        mesh.Cell1DsId.push_back(map1D_2[uM_key]);
                        mesh.Cell1DsExtrema.conservativeResize(2, map1D_2[uM_key] + 1);
                        mesh.Cell1DsExtrema(0, map1D_2[uM_key]) = c;
                        mesh.Cell1DsExtrema(1, map1D_2[uM_key]) = d;

                        archi.insert(uM_key);
                    }

                    unsigned int e = min(v,M_id);
                    unsigned int f = max(v,M_id);
                    pair<unsigned int,unsigned int> vM_key = {e,f};
                    if (!map1D_2.count(vM_key)) {
                        map1D_2[vM_key] = mesh.Cell1DsId.size();
                        mesh.Cell1DsId.push_back(map1D_2[vM_key]);
                        mesh.Cell1DsExtrema.conservativeResize(2, map1D_2[vM_key] + 1);
                        mesh.Cell1DsExtrema(0, map1D_2[vM_key]) = e;
                        mesh.Cell1DsExtrema(1, map1D_2[vM_key]) = f;

                        archi.insert(vM_key);
                    }
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
                    unsigned int a = min(bari_id[i], bari_id[j]);
                    unsigned int b = max(bari_id[i], bari_id[j]);
                    pair<unsigned int,unsigned int> BB_key = {a, b};
                    if (!map1D_2.count(BB_key)) {
                        map1D_2[BB_key] = mesh.Cell1DsId.size();
                        mesh.Cell1DsId.push_back(map1D_2[BB_key]);
                        mesh.Cell1DsExtrema.conservativeResize(2, map1D_2[BB_key] + 1);
                        mesh.Cell1DsExtrema(0, map1D_2[BB_key]) = a;
                        mesh.Cell1DsExtrema(1, map1D_2[BB_key]) = b;

                        archi.insert(BB_key);
                    }
                }
            }
        }
    }

    /// LATI NELLA MESH
    mesh.NumCell0Ds = mesh.Cell0DsId.size();   // NumCell0Ds
    mesh.NumCell1Ds = mesh.Cell1DsId.size();   // NumCell1Ds

    unordered_set<string> facce_codificate;

    // Ricostruzione mappa adiacenze
    unordered_map<unsigned int, unordered_set<unsigned int>> adiacenze;
    for (const auto& [u, v] : archi) {
        adiacenze[u].insert(v);
        adiacenze[v].insert(u);
    }

    for (const auto& [u, vicini_u] : adiacenze) {
        for (unsigned int v : vicini_u) {
            if (v <= u) continue;
            for (unsigned int w : adiacenze[v]) {
                if (w <= v || !adiacenze[u].count(w)) continue;
                if (u == v || v == w || w == u) continue;

                // Triangolo potenziale con vertici u,v,w
                // Determina ordine ciclico corretto: u → v → w → u
                vector<int> vertici = { (int)u, (int)v, (int)w };

                // Check se faccia già inserita (indipendentemente da ordine ciclico)
                vector<int> vertici_ordinati = vertici;
                sort(vertici_ordinati.begin(), vertici_ordinati.end());
                string codice = to_string(vertici_ordinati[0]) + "," + to_string(vertici_ordinati[1]) + "," + to_string(vertici_ordinati[2]);

                if (facce_codificate.count(codice)!=0) continue;
                    facce_codificate.insert(codice);

                // Costruzione degli archi nel giusto ordine
                vector<pair<int, int>> lati = {
                    {min(vertici[0], vertici[1]), max(vertici[0], vertici[1])},
                    {min(vertici[1], vertici[2]), max(vertici[1], vertici[2])},
                    {min(vertici[2], vertici[0]), max(vertici[2], vertici[0])}
                };

                vector<int> edges;
                bool arco_mancante = false;
                for (auto& lato : lati) {
                    if (map1D_2.count(lato)) {edges.push_back(map1D_2[lato]);}
                    else {
                        arco_mancante = true;
                        break;
                    }
                }
                if (arco_mancante) continue;

                mesh.Cell2DsId.push_back(fid);
                mesh.Cell2DsVertices.push_back(vertici);  // ordine ciclico coerente
                mesh.Cell2DsEdges.push_back(edges);       // ordine coerente con vertici
                mesh.Cell2DsNumEdges.push_back(3);
                fid++;
            }
        }
    }

    mesh.NumCell2Ds = mesh.Cell2DsId.size();   // NumCell2Ds

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


/*void build_duale(const PolyhedralMesh& geodetico, PolyhedralMesh& Goldberg) {

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
    vector<vector<int>> GoldbergFacce(geodetico.NumCell0Ds);   //per andare dal vertice del geodetico alla faccia di Goldberg

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
    

    // === Mappa: vertice -> facce che lo contengono ===
    map<int, vector<int>> verticeToFacce;
    for (unsigned int fid = 0; fid < geodetico.NumCell2Ds; ++fid) {
        for (int vid : geodetico.Cell2DsVertices[fid]) {
            verticeToFacce[vid].push_back(fid);
        }
    }

    // === Preprocess: arco orientato (v,u) -> faccia ===
    map<pair<int, int>, int> arcoToFaccia;
    for (unsigned int fid = 0; fid < geodetico.NumCell2Ds; ++fid) {
        const auto& vertici = geodetico.Cell2DsVertices[fid];
        int n = vertici.size();
        for (int i = 0; i < n; ++i) {
            int v = vertici[i];
            int u = vertici[(i + 1) % n];
            arcoToFaccia[{v, u}] = fid;
        }
    }

    // === Costruzione facce del duale inline ===
    vector<vector<int>> GoldbergFacce(geodetico.NumCell0Ds);

    for (unsigned int v = 0; v < geodetico.NumCell0Ds; ++v) {
        const auto& facce = verticeToFacce[v];
        if (facce.size() < 3) continue;

        // mappa: faccia -> vertici (ciclo)
        map<int, pair<int, int>> arcoPerFaccia; // (v,u) intorno a v
        for (int fid : facce) {
            const auto& verts = geodetico.Cell2DsVertices[fid];
            int n = verts.size();
            for (int i = 0; i < n; ++i) {
                if (verts[i] == v) {
                    int prev = verts[(i - 1 + n) % n];
                    int next = verts[(i + 1) % n];
                    arcoPerFaccia[fid] = {prev, next};
                    break;
                }
            }
        }

        // costruzione ciclo ordinato attorno al vertice v
        vector<int> ciclo;
        set<int> visitate;

        int fid_corrente = facce[0]; // iniziamo da una qualsiasi
        while (true) {
            if (visitate.count(fid_corrente)) break;
            visitate.insert(fid_corrente);
            ciclo.push_back(mapGoldbergVerticiID[fid_corrente]);

            int next_vert = arcoPerFaccia[fid_corrente].second;

            // cerca prossima faccia che ha (v, next_vert)
            bool trovato = false;
            for (int fid_prox : facce) {
                if (visitate.count(fid_prox)) continue;
                const auto& arco = arcoPerFaccia[fid_prox];
                if (arco.first == v && arco.second == next_vert) {
                    fid_corrente = fid_prox;
                    trovato = true;
                    break;
                }
            }

            if (!trovato) break;
        }

        if (ciclo.size() >= 3) {
            GoldbergFacce[v] = ciclo;
        }
    }


    // === Inserisci solo le facce valide nella mesh Goldberg ===
    Goldberg.Cell2DsId.clear();
    Goldberg.Cell2DsVertices.clear();
    Goldberg.Cell2DsNumEdges.clear();
    Goldberg.Cell2DsEdges.clear();

    for (const auto& face : GoldbergFacce) {
        if (face.size() >= 3) {
            Goldberg.Cell2DsId.push_back(Goldberg.Cell2DsId.size());
            Goldberg.Cell2DsVertices.push_back(face);
            Goldberg.Cell2DsNumEdges.push_back(face.size());
            Goldberg.Cell2DsEdges.push_back({});
        }
    }
    Goldberg.NumCell2Ds = Goldberg.Cell2DsId.size();


    /*Goldberg.NumCell2Ds = static_cast<unsigned int>(GoldbergFacce.size());
    Goldberg.Cell2DsId.resize(Goldberg.NumCell2Ds);
    Goldberg.Cell2DsVertices = GoldbergFacce;
    Goldberg.Cell2DsNumEdges.resize(Goldberg.NumCell2Ds);
    Goldberg.Cell2DsEdges.resize(Goldberg.NumCell2Ds);*/


    /*
    // === Cell1D: crea spigoli univoci dalle facce ===
    map<pair<int, int>, int> mapEdgeGoldberg;
    int eid = 0;
    for (unsigned int fid = 0; fid < Goldberg.NumCell2Ds; ++fid) {
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

    // === Costruzione spigoli ===
    map<pair<int, int>, int> mapEdgeGoldberg;
    int eid = 0;
        // Mappa: coppia (a,b) con orientamento -> ID spigolo
    map<pair<int, int>, int> orientedEdgeToId;
    eid = 0;

    for (unsigned int fid = 0; fid < Goldberg.NumCell2Ds; ++fid) {
        const auto& vertici = Goldberg.Cell2DsVertices[fid];
        size_t n = vertici.size();
        vector<int> edges_faccia;

        for (size_t i = 0; i < n; ++i) {
            int a = vertici[i];
            int b = vertici[(i + 1) % n];  // ciclo chiuso
            auto key = make_pair(a, b);   // orientato!

            int eid_corrente;
            if (orientedEdgeToId.count(key)) {
                eid_corrente = orientedEdgeToId[key];
            } else {
                eid_corrente = eid++;
                orientedEdgeToId[key] = eid_corrente;
                Goldberg.Cell1DsId.push_back(eid_corrente);
                Goldberg.Cell1DsExtrema.conservativeResize(2, eid);
                Goldberg.Cell1DsExtrema(0, eid_corrente) = a;
                Goldberg.Cell1DsExtrema(1, eid_corrente) = b;
            }

            edges_faccia.push_back(eid_corrente);
        }

        Goldberg.Cell2DsEdges[fid] = edges_faccia;
    }

    Goldberg.NumCell1Ds = Goldberg.Cell1DsId.size();


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
*/

void build_duale(const PolyhedralMesh& geodetico, PolyhedralMesh& duale) {
    duale = PolyhedralMesh();

    const int num_facce = geodetico.NumCell2Ds;

    // === 1. Costruisci i nuovi vertici: baricentri delle facce ===
    duale.NumCell0Ds = num_facce;
    duale.Cell0DsId.resize(num_facce);
    duale.Cell0DsCoordinates.resize(3, num_facce);

    for (int fid = 0; fid < num_facce; ++fid) {
        Vector3d baricentro(0, 0, 0);
        for (int vid : geodetico.Cell2DsVertices[fid]) {
            baricentro += geodetico.Cell0DsCoordinates.col(vid);
        }
        baricentro /= static_cast<double>(geodetico.Cell2DsVertices[fid].size());

        // Proiezione sulla sfera unitaria
        baricentro.normalize();

        duale.Cell0DsId[fid] = fid;
        duale.Cell0DsCoordinates.col(fid) = baricentro;
    }

    // === 2. Mappa vertice originale -> facce adiacenti ===
    unordered_map<int, vector<int>> vertice_to_facce;
    for (int fid = 0; fid < geodetico.NumCell2Ds; ++fid) {
        for (int vid : geodetico.Cell2DsVertices[fid]) {
            vertice_to_facce[vid].push_back(fid);
        }
    }

    // === 3. Costruisci le facce del duale (una per ogni vertice originale) ===
    const int num_facce_duale = geodetico.NumCell0Ds;
    duale.NumCell2Ds = num_facce_duale;
    duale.Cell2DsId.resize(num_facce_duale);
    duale.Cell2DsVertices.resize(num_facce_duale);
    duale.Cell2DsEdges.resize(num_facce_duale);
    duale.Cell2DsNumEdges.resize(num_facce_duale);


    for (int vid = 0; vid < num_facce_duale; ++vid) {
        const vector<int>& facce_adiacenti = vertice_to_facce[vid];

        // Ordinamento ciclico dei baricentri attorno al vertice
        vector<pair<double, int>> angoli;
        Vector3d center = geodetico.Cell0DsCoordinates.col(vid);

        for (int fid : facce_adiacenti) {
            Vector3d p = duale.Cell0DsCoordinates.col(fid);
            Vector3d v = p - center;
            double angle = atan2(v.y(), v.x());   //proiezione 2D semplificata
            angoli.emplace_back(angle, fid);
        }

        sort(angoli.begin(), angoli.end());

        vector<int> ordered_fids;
        for (auto& [angolo, fid] : angoli) {
            ordered_fids.push_back(fid);
        }

        duale.Cell2DsId[vid] = vid;
        duale.Cell2DsVertices[vid] = ordered_fids;
        duale.Cell2DsNumEdges[vid] = static_cast<int>(ordered_fids.size());
    }
    
    // === 3. Costruisci i lati del duale usando quelli del geodetico ===

    // Mappa da edge originale → facce adiacenti
    unordered_map<int, vector<int>> edge_to_faces;
    for (int fid = 0; fid < num_facce; ++fid) {
        for (int eid : geodetico.Cell2DsEdges[fid]) {
            edge_to_faces[eid].push_back(fid);
        }
    }

    // Per ogni edge condiviso da 2 facce, crea un lato nel duale
    map<pair<int, int>, int> edge_map;
    vector<pair<int, int>> edge_extrema;

    for (const auto& [eid, faces] : edge_to_faces) {
        if (faces.size() != 2) continue; // bordo o errore topologico
        int f1 = faces[0], f2 = faces[1];
        if (f1 == f2) continue;

        pair<int, int> key = minmax(f1, f2);
        if (edge_map.count(key)) continue;

        edge_map[key] = static_cast<int>(edge_extrema.size());
        edge_extrema.push_back(key);
    }

    duale.NumCell1Ds = static_cast<unsigned int>(edge_extrema.size());
    duale.Cell1DsId.resize(duale.NumCell1Ds);
    duale.Cell1DsExtrema = MatrixXi(2, duale.NumCell1Ds);

    for (int eid = 0; eid < duale.NumCell1Ds; ++eid) {
        duale.Cell1DsId[eid] = eid;
        duale.Cell1DsExtrema(0, eid) = edge_extrema[eid].first;
        duale.Cell1DsExtrema(1, eid) = edge_extrema[eid].second;
    }

    // === 4. Inserisci gli edge ID nelle facce del duale ===
    for (int fid = 0; fid < duale.NumCell2Ds; ++fid) {
        const vector<int>& v = duale.Cell2DsVertices[fid];
        vector<int> edge_ids;

        for (size_t i = 0; i < v.size(); ++i) {
            int from = v[i];
            int to = v[(i + 1) % v.size()];
            pair<int, int> key = minmax(from, to);
            auto it = edge_map.find(key);
            if (it != edge_map.end()) {
                edge_ids.push_back(it->second);
            }
        }

        duale.Cell2DsEdges[fid] = edge_ids;
    }

    // === 5. Definisci il singolo poliedro ===
    duale.NumCell3Ds = 1;
    duale.Cell3DsId = {0};

    vector<int> all_v(duale.NumCell0Ds), all_e(duale.NumCell1Ds), all_f(duale.NumCell2Ds);
    iota(all_v.begin(), all_v.end(), 0);
    iota(all_e.begin(), all_e.end(), 0);
    iota(all_f.begin(), all_f.end(), 0);

    duale.Cell3DsVertices = {all_v};
    duale.Cell3DsEdges    = {all_e};
    duale.Cell3DsFaces    = {all_f};
}

/*-----------------------------------------------------------------------------------------------*/

// Funzione per calcolare il cammino minimo

void Dijkstra(PolyhedralMesh& mesh, unsigned int vertice_iniziale, unsigned int vertice_finale)
{
    // Lista adiacenza dal vettore degli archi (Cell1Ds)
    map<unsigned int, vector<unsigned int>> ListaVerticiAdiacenti;   // key: id del vertice -> val: vettore dei vertici adiacenti in id

    for (unsigned int j = 0; j < mesh.NumCell1Ds; ++j) {
        unsigned int u = mesh.Cell1DsExtrema(0, j);
        unsigned int v = mesh.Cell1DsExtrema(1, j);
        ListaVerticiAdiacenti[u].push_back(v);
        ListaVerticiAdiacenti[v].push_back(u);
    }

    // Inizializza distanze e predecessori
    map<unsigned int, double> dist;   // key: id del vertice -> val: distanza dal vertice_iniziale
    map<unsigned int, unsigned int> pred;   // key: id del vertice -> val: id del suo vertice predecessore
    for (const auto& [vert, _] : ListaVerticiAdiacenti) {
        dist[vert] = numeric_limits<double>::infinity();   // inizialmente i vertici sono isolati rispetto al vertice_iniziale
    }
    dist[vertice_iniziale] = 0;   // il vertice iniziale ha distanza 0 rispetto a se stesso

    // coda di priorità: coppia (distanza, vertice), ordinata per distanza crescente
    priority_queue< pair<double, unsigned int>, vector<pair<double, unsigned int>>, greater<> > pq;
    pq.push({0.0, vertice_iniziale});

    // Algoritmo Dijkstra
    while (!pq.empty()) {
        double d = pq.top().first;
        unsigned int u = pq.top().second;
        pq.pop();

        if (d > dist[u]) continue;

        if (u == vertice_finale) break;  // se vuoi fermarti alla destinazione

        for (unsigned int v : ListaVerticiAdiacenti[u]) {
            // Calcola distanza euclidea tra u e v
            double dist_uv = sqrt(
                pow(mesh.Cell0DsCoordinates(0,u) - mesh.Cell0DsCoordinates(0,v), 2) +
                pow(mesh.Cell0DsCoordinates(1,u) - mesh.Cell0DsCoordinates(1,v), 2) +
                pow(mesh.Cell0DsCoordinates(2,u) - mesh.Cell0DsCoordinates(2,v), 2)
            );

            double nuova_dist = dist[u] + dist_uv;
            if (nuova_dist < dist[v]) {
                dist[v] = nuova_dist;
                pred[v] = u;
                pq.push({nuova_dist, v});
            }
        }
    }

    // Ricostruisci il cammino minimo
    vector<unsigned int> path;
    unsigned int current = vertice_finale;
    while (current != vertice_iniziale) {
        if (pred.find(current) == pred.end()) {
            cout << "Nessun cammino trovato!" << endl;
            return;
        }
        path.push_back(current);
        current = pred[current];
    }
    path.push_back(vertice_iniziale);
    reverse(path.begin(), path.end());

    // Salva risultati nella mesh
    mesh.num_archiPath = path.size() - 1;
    mesh.lunghezza_Path = dist[vertice_finale];

    // Marca i vertici nel cammino
    mesh.MarkerCell0Ds[1].clear();
    for (unsigned int v : path) {
        mesh.MarkerCell0Ds[1].push_back(v);
        cout << "Vertice cammino: " << v << endl;
    }

    // Marca gli archi nel cammino
    mesh.MarkerCell1Ds[1].clear();
    for (size_t i = 1; i < path.size(); ++i) {
        unsigned int u = path[i-1];
        unsigned int v = path[i];
        for (unsigned int j = 0; j < mesh.NumCell1Ds; ++j) {
            if ((mesh.Cell1DsExtrema(0, j) == u && mesh.Cell1DsExtrema(1, j) == v) ||
                (mesh.Cell1DsExtrema(0, j) == v && mesh.Cell1DsExtrema(1, j) == u)) {
                mesh.MarkerCell1Ds[1].push_back(mesh.Cell1DsId[j]);
                cout << "Arco cammino: " << j << endl;
                break;
            }
        }
    }
}

