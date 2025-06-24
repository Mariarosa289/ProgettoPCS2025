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
- function : normalizza
- input 1 : array con le 3 coordinate non normalizzate del punto p
- return : array con le 3 coordinate normalizzate del punto p
*/

array<double, 3> normalizza(const array<double, 3>& p) 
{
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);   // calcolo la norma del punto P
    return {p[0]/norm, p[1]/norm, p[2]/norm};   
}

/*-----------------------------------------------------------------------------------------------*/
/*
Funzione per arrotondare il punto, scegliendo un errore eps = 1e-6
- function : arrotonda_punto
- input 1 : array con le 3 coordinate del punto p
- input 2 : errore eps fissato a 1e-6
- return : array con le 3 coordinate arrotondate del punto p
*/

array<double, 3> arrotonda_punto(const array<double, 3>& p)
{
    double eps= 1e-6;
    return {
        round(p[0] / eps) * eps,
        round(p[1] / eps) * eps,
        round(p[2] / eps) * eps
    };
}

/*-----------------------------------------------------------------------------------------------*/
/*
Funzione che raccoglie i vertici e le facce iniziali di un icosaedro regolare centrato in origine e con i vertici a stessa distanza da esso
- function : build_ico
- input 1 : vector vuoto dei vertici in coordinate
- input 2 : vector vuoto delle facce in id dei loro vertici
- return : i vector di vertici e facce popolati
*/

void build_ico(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce) 
{
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
/*
Funzione che raccoglie i vertici e le facce iniziali di un ottaedro regolare centrato in origine e con i vertici a stessa distanza da esso
- function : build_octa
- input 1 : vector vuoto dei vertici in coordinate
- input 2 : vector vuoto delle facce in id dei loro vertici
- return : i vector di vertici e facce popolati
*/

void build_octa(vector<array<double, 3>>& vertici,vector<array<int, 3>>& facce) 
{
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
/*
Funzione che raccoglie i vertici e le facce iniziali di un tetraedro regolare centrato in origine e con i vertici a stessa distanza da esso
- function : build_tetra
- input 1 : vector vuoto dei vertici in coordinate
- input 2 : vector vuoto delle facce in id dei loro vertici
- return : i vector di vertici e facce popolati
*/

void build_tetra(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce) 
{
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
/*
Funzione che memorizza il vertice normalizzato
- function : salva_vertice_norm
- input 1 : array con le coordinate del vertice p 
- input 2 : mappa map0D in cui salvare le coordinate normalizzate con l'id del vertice
- input 3 : mesh da popolare
- return : vid aumentato di uno 
*/

unsigned int salva_vertice_norm(const array<double,3>& coordinate,
                                map<array<double,3>, unsigned int>& map0D,
                                PolyhedralMesh& mesh) 
{
    static unsigned int vid = 0;   // static: per poter aumentare il vid ad ogni chiamata

    auto vertice_norma = normalizza(coordinate);   // normalizza vertice
    auto vertice_arrotondato = arrotonda_punto(vertice_norma);   // arrotonda vertice

    // vediamo se il vertice è già nella map0D
    auto it = map0D.find(vertice_arrotondato);
    if (it != map0D.end()) return it -> second;  

    // vertice mancante nella map0D => aggiungere nella map0D e nella mesh
    map0D[vertice_arrotondato] = vid;  
    mesh.Cell0D_id.push_back(vid);
    mesh.Cell0D_coordinate.conservativeResize(3, vid + 1);   // conservativeResize : resize mantendendo i dati già presenti
    for (int i = 0; i < 3; ++i)
        mesh.Cell0D_coordinate(i, vid) = vertice_norma[i];

    return vid++;
}


/*-----------------------------------------------------------------------------------------------*/
/*
Funzione che popola il mesh dopo aver applicato la triangolazione 1
- function : build_classe_1
- input 1 : p : numero di vertici che si osserva guardando ciascuna faccia del poligono regolare iniziale
- input 2 : q : valenza (numero di facce adiacenti a un vertice del poligono regolare iniziale)
- input 3 : b : numero di divisione di ogni lato di una faccia base
- input 4 : c : classe
- input 5 : mesh da popolare
- return : popolamento del mesh 
*/

void build_classe_1(unsigned int p, 
                    unsigned int q, 
                    unsigned int b, 
                    unsigned int c, 
                    PolyhedralMesh& mesh) 
{
    /// Vertici e facce del poligono regolare di base
    vector<array<double, 3>> vertici_iniziali;
    vector<array<int, 3>> facce_iniziali;

    switch (q) {
        case 5:
            build_ico(vertici_iniziali, facce_iniziali);
            break;
        case 4:
            build_octa(vertici_iniziali, facce_iniziali);
            break;
        case 3:
            build_tetra(vertici_iniziali, facce_iniziali);
            break;
    }

    map<array<double, 3>, unsigned int> map0D;   // raccolta univoca dei vertici. key = coordinate del vertice -> id del vertice
    map<pair<unsigned int, unsigned int>, unsigned int> map1D;   // raccolta univoca dei lati. key = id dei vertici estremi -> id del lato
    unsigned int fid = 0;   // indice di faccia

    /// Triangolazione 1 per ogni faccia base + popolamento del mesh
    for (const auto& faccia : facce_iniziali) {
        const auto& A = vertici_iniziali[faccia[0]];
        const auto& B = vertici_iniziali[faccia[1]];
        const auto& C = vertici_iniziali[faccia[2]];

        /// popolamento mesh di vertici della triangolazione
        vector<vector<unsigned int>> griglia(b + 1);   // griglia triangolare
        for (unsigned int i = 0; i <= b; ++i) {
            griglia[i].resize(i + 1);   
            for (unsigned int j = 0; j <= i; ++j) {   // ATTENZIONE: FORMULA?
                double u = 1.0 - static_cast<double>(i) / b;   // static_cast per conversione sicura
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

        /// popolamento del resto del mesh
        for (unsigned int i = 0; i < b; ++i) {
            for (unsigned int j = 0; j < i + 1; ++j) {
                unsigned int v1 = griglia[i][j];
                unsigned int v2 = griglia[i+1][j];
                unsigned int v3 = griglia[i+1][j + 1];

                mesh.Cell2D_id.push_back(fid++);   // popolamento Cell2D_id
                mesh.Cell2D_vertici.push_back({(int)v1, (int)v2, (int)v3});   // popolamento Cell2D_vertici

                vector<int> lati_diFaccia;   // raccoglie i lati della faccia

                for (int k = 0; k < 3; ++k) {
                    unsigned int a = mesh.Cell2D_vertici.back()[k];
                    unsigned int b = mesh.Cell2D_vertici.back()[(k + 1) % 3];
                    auto lato = make_pair(a, b);
                    auto lato_inverso = make_pair(b, a);
                    
                    if (map1D.count(lato) == 0) {
                        if (map1D.count(lato_inverso) == 0) {   // lato e lato_inverso non presenti -> ce lo aggiungo
                            map1D[lato] = mesh.Cell1D_id.size();
                            mesh.Cell1D_id.push_back(map1D[lato]);   // popolamento Cell1D_id
                            mesh.Cell1D_estremi.conservativeResize(2, map1D[lato] + 1);   // popolamento Cell1D_estremi
                            mesh.Cell1D_estremi(0, map1D[lato]) = lato.first;
                            mesh.Cell1D_estremi(1, map1D[lato]) = lato.second; 
                            lati_diFaccia.push_back(map1D[lato]);
                        } else if (map1D.count(lato_inverso) != 0) {
                            lati_diFaccia.push_back(map1D[lato_inverso]);
                        }
                    } else if(map1D.count(lato) != 0) {
                        lati_diFaccia.push_back(map1D[lato]);
                    }
                }

                mesh.Cell2D_lati.push_back(lati_diFaccia);   // popolamento Cell2D_lati
                mesh.Cell2D_numLati.push_back(lati_diFaccia.size());   //popolamento Cell2D_numLati

                /// stesso popolamento, ma nel caso j < 1
                if (j < i) {   //ATTENZIONE: che dice la condizione?
                    unsigned int v4 = griglia[i][j+1];
                    mesh.Cell2D_id.push_back(fid++);  
                    mesh.Cell2D_vertici.push_back({(int)v1, (int)v3, (int)v4});  

                    vector<int> lati_diFaccia_2;

                    for (int k = 0; k < 3; ++k) {
                    unsigned int a = mesh.Cell2D_vertici.back()[k];
                    unsigned int b = mesh.Cell2D_vertici.back()[(k + 1) % 3];
                    auto lato = make_pair(a, b);
                    auto lato_inverso = make_pair(b, a);
                    
                    if (map1D.count(lato) == 0) {
                        if (map1D.count(lato_inverso) == 0) {   
                            map1D[lato] = mesh.Cell1D_id.size();
                            mesh.Cell1D_id.push_back(map1D[lato]);
                            mesh.Cell1D_estremi.conservativeResize(2, map1D[lato] + 1);
                            mesh.Cell1D_estremi(0, map1D[lato]) = lato.first;
                            mesh.Cell1D_estremi(1, map1D[lato]) = lato.second;
                            lati_diFaccia_2.push_back(map1D[lato]);
                        } else if (map1D.count(lato_inverso) != 0) {
                            lati_diFaccia_2.push_back(map1D[lato_inverso]);
                        }
                    } else if(map1D.count(lato) != 0) {
                        lati_diFaccia_2.push_back(map1D[lato]);
                    }
                }

                mesh.Cell2D_lati.push_back(lati_diFaccia_2);
                mesh.Cell2D_numLati.push_back(3);
                }
            }
        }
    }

    mesh.Cell0D_num = mesh.Cell0D_id.size();   // popolamento Cell0D_num
    mesh.Cell1D_num = mesh.Cell1D_id.size();   // popolamento Cell1D_num
    mesh.Cell2D_num = mesh.Cell2D_id.size();   // popolamento Cell2D_num

    /// Popolamento della Cell3D
    mesh.Cell3D_num = 1;
    mesh.Cell3D_id.push_back(0);

    vector<int> vertici_cell3d(mesh.Cell0D_num);
    iota(vertici_cell3d.begin(), vertici_cell3d.end(), 0);
    mesh.Cell3D_vertici.push_back(vertici_cell3d);

    vector<int> lati_cell3d(mesh.Cell1D_num);
    iota(lati_cell3d.begin(), lati_cell3d.end(), 0);
    mesh.Cell3D_lati.push_back(lati_cell3d);

    vector<int> facce_cell3d(mesh.Cell2D_num);
    iota(facce_cell3d.begin(), facce_cell3d.end(), 0);
    mesh.Cell3D_facce.push_back(facce_cell3d);

    // Ottimizzazione memoria : togliere memoria inutilizzata dopo il dinamismo
    mesh.Cell2D_id.shrink_to_fit();
    mesh.Cell2D_lati.shrink_to_fit();
    mesh.Cell2D_vertici.shrink_to_fit();
    mesh.Cell2D_numLati.shrink_to_fit();  

    // ATTENZIONE: fare lo shrink_to_fit anche per le altre Cell?
}


/*-----------------------------------------------------------------------------------------------*/
/*
- Funzione : baricentro : calcola in coordinate il baricentro di un triangolo
- Inputs : A,B,C : vertici del triangolo (in coord)
- return : coordinate del baricentro di quel triangolo
*/
array<double,3> calcola_baricentro(const array<double,3>& A, 
                                   const array<double,3>& B, 
                                   const array<double,3>& C) 
{ return {(A[0]+B[0]+C[0])/3.0, (A[1]+B[1]+C[1])/3.0, (A[2]+B[2]+C[2])/3.0}; }


/*-----------------------------------------------------------------------------------------------*/
/*
- Funzione : calcola_puntomedio : calcola in coordinate  il punto medio di un lato
- Inputs : A,B : estremi di un lato (in coord)
- return : punto medio in coordinate
*/
array<double,3> calcola_puntomedio(const array<double,3>& A, const array<double,3>& B) 
{ return {(A[0]+B[0])/2.0, (A[1]+B[1])/2.0, (A[2]+B[2])/2.0}; }


/*-----------------------------------------------------------------------------------------------*/
/* 
- Funzione : build_classe_2 : mesh secondo la triangolazione di classe 2 con b = c con b>0 (classe II)
- Inputs : p, q, b, c==b con b > 0
- return : riempimento del mesh in PolyhedralMesh
*/
void build_classe_2(unsigned int p, 
                    unsigned int q, 
                    unsigned int b, 
                    unsigned int c, 
                    PolyhedralMesh& mesh) {
    
    /// vertici e facce del poligono base
    vector<array<double, 3>> vertici_iniziali;
    vector<array<int, 3>> facce_iniziali;

    switch (q) {
        case 5:
            build_ico(vertici_iniziali, facce_iniziali);
            break;
        case 4:
            build_octa(vertici_iniziali, facce_iniziali);
            break;
        case 3:
            build_tetra(vertici_iniziali, facce_iniziali);
            break;
    }

    map<array<double, 3>, unsigned int> map0D;   // key = coordinate del vertice -> ide del vertice
    map<pair<unsigned int,unsigned int>, unsigned int> map1D;   // key = id dei vertici estremi -> id del lato
    
    unsigned int fid = 0;   // id delle facce

    set<pair<unsigned int,unsigned int>> archi;   // raccolta di tutti gli archi identificati dagli id dei vertici

    for (const auto& faccia : facce_iniziali) {
        const auto& A = vertici_iniziali[faccia[0]];
        const auto& B = vertici_iniziali[faccia[1]];
        const auto& C = vertici_iniziali[faccia[2]]; 

        map<array<int,3>, unsigned int> coordinate_bari_id;   // key: array con coord bari -> id del bari  
        coordinate_bari_id.clear();

        /// "griglia" come nel build_classe_1
        for (int i = 0; i <= (int)b; ++i) {
            for (int j = 0; j <= (int)(b - i); ++j) {
                int k = b - i - j;   // regola base delle coord baricentrice
                double x = (i * A[0] + j * B[0] + k * C[0]) / double(b);   //combinaz convessa del triangolo iniziale
                double y = (i * A[1] + j * B[1] + k * C[1]) / double(b);
                double z = (i * A[2] + j * B[2] + k * C[2]) / double(b);
                array<double,3> P = {x,y,z};
                coordinate_bari_id[{i,j,k}] = salva_vertice_norm(P, map0D, mesh);  
            }
        }


        vector<tuple<unsigned int,unsigned int,unsigned int>> sotto_triangoli;   // vettore di sottotriangoli rappresentati dagli id

        /// Salvataggio di Sotto-triangoli
        for (int i = 0; i < (int)b; ++i) {
            for (int j = 0; j < (int)(b - i); ++j) {
                int k = b - i - j;
                unsigned int v1 = coordinate_bari_id[{i,j,k}];
                unsigned int v2 = coordinate_bari_id[{i+1,j,k-1}];
                unsigned int v3 = coordinate_bari_id[{i,j+1,k-1}];
                sotto_triangoli.emplace_back(v1, v2, v3);   // emplace_back : costruisce l'elemento direttamente dentro il vector

                if (i + j < (int)(b - 1)) {
                    unsigned int v4 = coordinate_bari_id[{i+1,j+1,k-2}];
                    sotto_triangoli.emplace_back(v2, v4, v3);
                }
            }
        }


        /// Lati condivisi da due sotto-triangoli
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

        vector<unsigned int> bari_id;   // vector di id dei baricentri

        /// Analisi di ogni sotto-triangolo
        for (auto& tri : sotto_triangoli) {   

            unsigned int v1 = get<0>(tri); 
            unsigned int v2 = get<1>(tri);
            unsigned int v3 = get<2>(tri);

            /// Calcolo del baricentro di ogni sotto-triangolo
            array<double,3> A = {mesh.Cell0D_coordinate(0,v1), mesh.Cell0D_coordinate(1,v1), mesh.Cell0D_coordinate(2,v1)};
            array<double,3> B = {mesh.Cell0D_coordinate(0,v2), mesh.Cell0D_coordinate(1,v2), mesh.Cell0D_coordinate(2,v2)};
            array<double,3> C = {mesh.Cell0D_coordinate(0,v3), mesh.Cell0D_coordinate(1,v3), mesh.Cell0D_coordinate(2,v3)};
            
            array<double,3> G = calcola_baricentro(A,B,C);   // calcolo bari del sottotriangolo
            unsigned int G_id = salva_vertice_norm(G, map0D, mesh);   // Cell0D_coordinate e Cell0D_id del baricentro

            bari_id.push_back(G_id);

            /// lati bari-vertici_sottotriangolo
            for (auto v : {v1, v2, v3}) {
                unsigned int a = min(G_id,v);
                unsigned int b = max(G_id,v);
                pair<unsigned int,unsigned int> GV_key = {a,b};
                
                if (!map1D.count(GV_key)) {
                    map1D[GV_key] = mesh.Cell1D_id.size();
                    mesh.Cell1D_id.push_back(map1D[GV_key]);
                    mesh.Cell1D_estremi.conservativeResize(2, map1D[GV_key] + 1);
                    mesh.Cell1D_estremi(0, map1D[GV_key]) = a;
                    mesh.Cell1D_estremi(1, map1D[GV_key]) = b;

                    archi.insert(GV_key);
                }
            }

            /// lati coi punti medi
            array<pair<unsigned int,unsigned int>,3> lati = {make_pair(v1,v2),make_pair(v2,v3),make_pair(v3,v1)};   //lati del sottotriangolo
        
            for (auto [u,v] : lati) {
                if (u > v) swap(u,v);   //  ordine (min, max)

                /// punti medi
                if (lati_condivisi[{u,v}] == 1) {
                    array<double,3> P1 = {mesh.Cell0D_coordinate(0,u), mesh.Cell0D_coordinate(1,u), mesh.Cell0D_coordinate(2,u)};
                    array<double,3> P2 = {mesh.Cell0D_coordinate(0,v), mesh.Cell0D_coordinate(1,v), mesh.Cell0D_coordinate(2,v)};
                    array<double,3> M = calcola_puntomedio(P1, P2);
                    unsigned int M_id = salva_vertice_norm(M, map0D, mesh);

                    /// lati bari-punto_medio
                    unsigned int a = min(M_id,G_id);
                    unsigned int b = max(M_id,G_id);
                    pair<unsigned int,unsigned int> MG_key = {a,b};

                    if (!map1D.count(MG_key)) {
                        map1D[MG_key] = mesh.Cell1D_id.size();
                        mesh.Cell1D_id.push_back(map1D[MG_key]);
                        mesh.Cell1D_estremi.conservativeResize(2, map1D[MG_key] + 1);
                        mesh.Cell1D_estremi(0, map1D[MG_key]) = a;
                        mesh.Cell1D_estremi(1, map1D[MG_key]) = b;

                        archi.insert(MG_key);
                    }

                    /// lati punto_medio-vertici_generatori
                    unsigned int c = min(u,M_id);
                    unsigned int d = max(u,M_id);
                    pair<unsigned int,unsigned int> uM_key = {c,d};

                    if (!map1D.count(uM_key)) {
                        map1D[uM_key] = mesh.Cell1D_id.size();
                        mesh.Cell1D_id.push_back(map1D[uM_key]);
                        mesh.Cell1D_estremi.conservativeResize(2, map1D[uM_key] + 1);
                        mesh.Cell1D_estremi(0, map1D[uM_key]) = c;
                        mesh.Cell1D_estremi(1, map1D[uM_key]) = d;

                        archi.insert(uM_key);
                    }

                    unsigned int e = min(v,M_id);
                    unsigned int f = max(v,M_id);
                    pair<unsigned int,unsigned int> vM_key = {e,f};

                    if (!map1D.count(vM_key)) {
                        map1D[vM_key] = mesh.Cell1D_id.size();
                        mesh.Cell1D_id.push_back(map1D[vM_key]);
                        mesh.Cell1D_estremi.conservativeResize(2, map1D[vM_key] + 1);
                        mesh.Cell1D_estremi(0, map1D[vM_key]) = e;
                        mesh.Cell1D_estremi(1, map1D[vM_key]) = f;

                        archi.insert(vM_key);
                    }
                }
            }
        } 

        /// lati tra di sottotriangoli adiacenti
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

                    if (!map1D.count(BB_key)) {
                        map1D[BB_key] = mesh.Cell1D_id.size();
                        mesh.Cell1D_id.push_back(map1D[BB_key]);
                        mesh.Cell1D_estremi.conservativeResize(2, map1D[BB_key] + 1);
                        mesh.Cell1D_estremi(0, map1D[BB_key]) = a;
                        mesh.Cell1D_estremi(1, map1D[BB_key]) = b;

                        archi.insert(BB_key);
                    }
                }
            }
        }
    }

    mesh.Cell0D_num = mesh.Cell0D_id.size();   // popolamento Cell0D_num
    mesh.Cell1D_num = mesh.Cell1D_id.size();   // popolamento Cell1D_num

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
                    if (map1D.count(lato)) {edges.push_back(map1D[lato]);}
                    else {
                        arco_mancante = true;
                        break;
                    }
                }
                if (arco_mancante) continue;

                mesh.Cell2D_id.push_back(fid);   //popolamento Cell2D_id
                mesh.Cell2D_vertici.push_back(vertici);  // popolamento Cell2D_vertici
                mesh.Cell2D_lati.push_back(edges);       // popolamento Cell2D_lati
                mesh.Cell2D_numLati.push_back(edges.size());   // popolamento Cell2D_numLati
                fid++;
            }
        }
    }

    mesh.Cell2D_num = mesh.Cell2D_id.size();   // popolamento Cell2D_num

	/// popolamento Cell3D
    mesh.Cell3D_num = 1;
    mesh.Cell3D_id.push_back(0);

    // Tutti i vertici usati nella mesh
    vector<int> vertici_cell3d(mesh.Cell0D_num);
    iota(vertici_cell3d.begin(), vertici_cell3d.end(), 0);
    mesh.Cell3D_vertici.push_back(vertici_cell3d);

    // Tutti gli spigoli usati nella mesh
    vector<int> lati_cell3d(mesh.Cell1D_num);
    iota(lati_cell3d.begin(), lati_cell3d.end(), 0);
    mesh.Cell3D_lati.push_back(lati_cell3d);

    // Tutte le facce usate nella mesh
    vector<int> facce_cell3d(mesh.Cell2D_num);
    iota(facce_cell3d.begin(), facce_cell3d.end(), 0);
    mesh.Cell3D_facce.push_back(facce_cell3d);

    /// ottimizzazione della memoria
    mesh.Cell2D_id.shrink_to_fit();
    mesh.Cell2D_lati.shrink_to_fit();
    mesh.Cell2D_vertici.shrink_to_fit();
    mesh.Cell2D_numLati.shrink_to_fit();

    // ATTENZIONE : valutare anche il shrink_to_fit per le altre Cell?
}


/*-----------------------------------------------------------------------------------------------*/
/*
Funzione che genera tutti i file .txt
*/
void genera_tutti_file(const PolyhedralMesh& mesh, 
                     const string& outfilename0D,
                     const string& outfilename1D,
                     const string& outfilename2D,
                     const string& outfilename3D) 
{
    if (! genera_file_Cell0D (mesh, outfilename0D))
    { 
        cerr<<"Problemi di esportazione per Cell0Ds.txt"<< endl;
    } else 
        cout<<"Esportazione terminata con successo per Cell0Ds.txt"<< endl;
	 
	if (! genera_file_Cell1D (mesh, outfilename1D))
	{
        cerr<<"Problemi di esportazione per Cell1Ds.txt"<< endl;
	} else
        cout<<"Esportazione terminata con successo per Cell1Ds.txt"<< endl;
	  	
	if (! genera_file_Cell2D (mesh, outfilename2D))
	{
        cerr<<"Problemi di esportazione per Cell2Ds.txt"<< endl;
	} else
        cout<<"Esportazione terminata con successo per Cell2Ds.txt"<< endl;
	  	
	  	
	if (! genera_file_Cell3D (mesh, outfilename3D))
	{
        cerr<<"Problemi di esportazione per Cell3Ds.txt"<< endl;
	} else
        cout<<"Esportazione terminata con successo per Cell3Ds.txt"<< endl;
	  	
}
		 
/*-----------------------------------------------------------------------------------------------*/
/*
Funzione che genera il file Cell0D
*/

bool genera_file_Cell0D(const PolyhedralMesh& mesh, const string& outfilename) {
    ofstream file(outfilename);
    if (!file.is_open()) return false;

    file << "Id, x, y, z\n";
    
    // controllo marker (non obbligatorio)
    auto it = mesh.Cell0D_marker.find(1);
    if (it != mesh.Cell0D_marker.end()) {
        cout << "La chiave 1 è presente nella mappa." << endl;
    } else {
        cout << "La chiave 1 non è presente nella mappa." << endl;
    }


    for (unsigned int i = 0; i < mesh.Cell0D_num; ++i) {   //stampa sul outputFile
        file << mesh.Cell0D_id[i] << ", "
             << setprecision(10) << mesh.Cell0D_coordinate(0, i) << ", "
             << setprecision(10) << mesh.Cell0D_coordinate(1, i) << ", "
			 << setprecision(10) << mesh.Cell0D_coordinate(2, i) << "\n";
            
    }

    file.close();
    
    return true;
}

/*-----------------------------------------------------------------------------------------------*/
/*
Funzione che genera il file Cell1D
*/

bool genera_file_Cell1D(const PolyhedralMesh& mesh, const string& outfilename) {
    ofstream file(outfilename);
    if (!file.is_open()) {
        cerr << "Errore apertura file " << outfilename << endl;
        return false;
    }

    file << "Id, Id_Origine, Id_Fine\n";

    // controllo marker (non obbligatorio)
    auto it = mesh.Cell1D_marker.find(1);
    if (it != mesh.Cell1D_marker.end()) {
        cout << "La chiave 1 è presente nella mappa." << endl;
    } else {
        cout << "La chiave 1 non è presente nella mappa." << endl;
    }

    for (unsigned int i= 0; i < mesh.Cell1D_num; ++i) {   // stampa sul outputFile
        file << mesh.Cell1D_id[i] << ", " 
		<< setprecision(10) << mesh.Cell1D_estremi(0,i)<< ", " 
		<<setprecision(10) << mesh.Cell1D_estremi(1,i) << "\n";
    }

    file.close();

    return true;
	
}

/*-----------------------------------------------------------------------------------------------*/
/*
Funzione che genera il file Cell2D
*/

bool genera_file_Cell2D(const PolyhedralMesh& mesh, const string& outfilename) {
    ofstream file(outfilename);
    if (!file.is_open()) {
        cerr << "Errore apertura file " << outfilename << endl;
        return false;
    }

    file << "Id, NumVertici, NumSpigoli, Id_Vertici, Id_Spigoli\n";

    for (unsigned int i= 0; i < mesh.Cell2D_num; ++i) {
        file << mesh.Cell2D_id[i] << ", "
		<< setprecision(10) << mesh.Cell2D_vertici[i].size() << ", "
        << setprecision(10) << mesh.Cell2D_lati[i].size() << ", ";
		
		file << "[" ;
        // Vertici
        for (size_t j=0; j < mesh.Cell2D_vertici[i].size();j++)
            file << mesh.Cell2D_vertici[i][j] << "  ";
        // Spigoli
        file << "], [" ;
        for (size_t j=0; j<mesh.Cell2D_lati[i].size();j++)
            file << mesh.Cell2D_lati[i][j] << "  ";

        file << "]\n";
    }

    file.close();

    return true;
}

/*-----------------------------------------------------------------------------------------------*/
/*
Funzione che genera Cell3D
*/
bool genera_file_Cell3D(const PolyhedralMesh& mesh, const string& outfilename) {
    ofstream file(outfilename);
    if (!file.is_open()) {
        cerr << "Errore apertura file " << outfilename << endl;
        return false;
    }

    file << "Id, NumVertici, NumSpigoli, NumFacce, Id_Vertici, Id_Spigoli, Id_Facce\n";
	
	for (unsigned int i = 0; i < mesh.Cell3D_num; ++i) {
        file << mesh.Cell0D_id[i] << ", "
		<< setprecision(10) << mesh.Cell3D_vertici[i].size() << ", "
		<< setprecision(10) << mesh.Cell3D_lati[i].size() << ", "
		<< setprecision(10) << mesh.Cell3D_facce[i].size() << ", ";
		
		file << "[" ;
        // Vertici
        for (size_t j=0; j<mesh.Cell3D_vertici[i].size();j++)
            file << mesh.Cell3D_vertici[i][j] << "  ";
		
		// Spigoli
        file << "], [" ;
        for (size_t j=0; j<mesh.Cell3D_lati[i].size();j++)
            file << mesh.Cell3D_lati[i][j] << "  ";
		
		// Facce
        file << "], [" ;
        for (size_t j=0; j<mesh.Cell3D_facce[i].size();j++)
            file << mesh.Cell3D_facce[i][j] << "  ";

        file << "]\n";
		
    }

    file.close();

    return true;
}


/*-----------------------------------------------------------------------------------------------*/
/*
Function: calcola il calcola_coordinate_centroide
input2: Cell0D_coordinate
*/
static Vector3d calcola_coordinate_centroide(const vector<int>& VerticesId, const MatrixXd& Cell0D_coordinate) {
    Vector3d c(0.0, 0.0, 0.0);

    for (int vid : VerticesId) {   // itera su ogni punto (con coordinate in colonna)
        c += Cell0D_coordinate.col(vid);
    }
    c /= static_cast<double>(VerticesId.size());   // media (statica_cast trasforma l'intero in double per evitare la diviosione intera)
    
	double norm = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);   // calcolo la norma per normalizzare
    Vector3d p(c[0]/norm, c[1]/norm, c[2]/norm);
	
    return p;     
}


/*-----------------------------------------------------------------------------------------------*/
/* Funzione che genera il duale partendo da un geodetico
- Funzione : build_duale
- input 1 : mesh geodetico da cui partire
- input 2 : mesh del duale
*/
void build_duale(const PolyhedralMesh& geodetico, PolyhedralMesh& duale) {
    duale = PolyhedralMesh();

    const int num_facce = geodetico.Cell2D_num;

    /// Costruisci i vertici del duale: baricentri delle facce del geodetico
    duale.Cell0D_num = num_facce;
    duale.Cell0D_id.resize(num_facce);
    duale.Cell0D_coordinate.resize(3, num_facce);

    for (int fid = 0; fid < num_facce; ++fid) {
        Vector3d baricentro(0, 0, 0);
        for (int vid : geodetico.Cell2D_vertici[fid]) {
            baricentro += geodetico.Cell0D_coordinate.col(vid);
        }
        baricentro /= static_cast<double>(geodetico.Cell2D_vertici[fid].size());

        // Proiezione sulla sfera unitaria
        baricentro.normalized();

        duale.Cell0D_id[fid] = fid;
        duale.Cell0D_coordinate.col(fid) = baricentro;
    }

    /// Mappa vertice originale -> facce adiacenti 
    unordered_map<int, vector<int>> vertice_to_facce;

    for (int fid = 0; fid < geodetico.Cell2D_num; ++fid) {
        for (int vid : geodetico.Cell2D_vertici[fid]) {
            vertice_to_facce[vid].push_back(fid);
        }
    }

    /// Costruisci le facce del duale (una per ogni vertice originale)
    const int num_facce_duale = geodetico.Cell0D_num;
    duale.Cell2D_num = num_facce_duale;
    duale.Cell2D_id.resize(num_facce_duale);
    duale.Cell2D_vertici.resize(num_facce_duale);
    duale.Cell2D_lati.resize(num_facce_duale);
    duale.Cell2D_numLati.resize(num_facce_duale);


    for (int vid = 0; vid < num_facce_duale; ++vid) {
        const vector<int>& facce_adiacenti = vertice_to_facce[vid];

        // Ordinamento ciclico dei baricentri attorno al vertice
        vector<pair<double, int>> angoli;
        Vector3d center = geodetico.Cell0D_coordinate.col(vid);

        for (int fid : facce_adiacenti) {
            Vector3d p = duale.Cell0D_coordinate.col(fid);
            Vector3d v = p - center;
            double angle = atan2(v.y(), v.x());   //proiezione 2D semplificata
            angoli.emplace_back(angle, fid);
        }

        sort(angoli.begin(), angoli.end());

        vector<int>fid_ordinate;
        for (auto& [angolo, fid] : angoli) {
           fid_ordinate.push_back(fid);
        }

        duale.Cell2D_id[vid] = vid;
        duale.Cell2D_vertici[vid] = fid_ordinate;
        duale.Cell2D_numLati[vid] = static_cast<int>(fid_ordinate.size());
    }
    
    /// Costruiamo i lati del duale usando quelli del geodetico
    // Mappa da edge originale → facce adiacenti
    unordered_map<int, vector<int>> lati_to_facce;

    for (int fid = 0; fid < num_facce; ++fid) {
        for (int eid : geodetico.Cell2D_lati[fid]) {
            lati_to_facce[eid].push_back(fid);
        }
    }

    // Per ogni lato condiviso da 2 facce, crea un lato nel duale
    map<pair<int, int>, int> lati_map;
    vector<pair<int, int>> estremi_lati;

    for (const auto& [eid, faces] : lati_to_facce) {
        if (faces.size() != 2) continue;   // bordo o errore topologico
        int f1 = faces[0], f2 = faces[1];
        if (f1 == f2) continue;

        pair<int, int> key = minmax(f1, f2);
        if (lati_map.count(key)) continue;

        lati_map[key] = static_cast<int>(estremi_lati.size());
        estremi_lati.push_back(key);
    }

    duale.Cell1D_num = static_cast<unsigned int>(estremi_lati.size());
    duale.Cell1D_id.resize(duale.Cell1D_num);
    duale.Cell1D_estremi = MatrixXi(2, duale.Cell1D_num);

    for (int eid = 0; eid < duale.Cell1D_num; ++eid) {
        duale.Cell1D_id[eid] = eid;
        duale.Cell1D_estremi(0, eid) = estremi_lati[eid].first;
        duale.Cell1D_estremi(1, eid) = estremi_lati[eid].second;
    }

    /// Inseriamo gli id dei lati nelle facce
    for (int fid = 0; fid < duale.Cell2D_num; ++fid) {
        const vector<int>& v = duale.Cell2D_vertici[fid];
        vector<int> lati_id;

        for (size_t i = 0; i < v.size(); ++i) {
            int from = v[i];
            int to = v[(i + 1) % v.size()];
            pair<int, int> key = minmax(from, to);
            auto it = lati_map.find(key);
            if (it != lati_map.end()) {
                lati_id.push_back(it->second);
            }
        }
        duale.Cell2D_lati[fid] = lati_id;
    }

    /// popolamento Cell3D
    duale.Cell3D_num = 1;
    duale.Cell3D_id = {0};

    vector<int> tutti_vertici(duale.Cell0D_num), tutti_lati(duale.Cell1D_num), tutte_facce(duale.Cell2D_num);
    iota(tutti_vertici.begin(), tutti_vertici.end(), 0);
    iota(tutti_lati.begin(), tutti_lati.end(), 0);
    iota(tutte_facce.begin(), tutte_facce.end(), 0);

    duale.Cell3D_vertici = {tutti_vertici};
    duale.Cell3D_lati    = {tutti_lati};
    duale.Cell3D_facce    = {tutte_facce};
}

/*-----------------------------------------------------------------------------------------------*/
/* 
- funzione : Dijkstra
*/

void Dijkstra(PolyhedralMesh& mesh, unsigned int vertice_iniziale, unsigned int vertice_finale)
{
    // Lista adiacenza dal vettore degli archi (Cell1Ds)
    map<unsigned int, vector<unsigned int>> ListaVerticiAdiacenti;   // key: id del vertice -> val: vettore dei vertici adiacenti in id

    for (unsigned int j = 0; j < mesh.Cell1D_num; ++j) {
        unsigned int u = mesh.Cell1D_estremi(0, j);
        unsigned int v = mesh.Cell1D_estremi(1, j);
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

        if (d > dist[u]) continue;   // controllo se la distanza analizzata da pq è peggiore della distanza salvata in dist

        if (u == vertice_finale) break;   // ci fermiama al vertice finale

        for (unsigned int v : ListaVerticiAdiacenti[u]) {
            // Calcola distanza euclidea tra u e v
            double dist_uv = sqrt(
                pow(mesh.Cell0D_coordinate(0,u) - mesh.Cell0D_coordinate(0,v), 2) +
                pow(mesh.Cell0D_coordinate(1,u) - mesh.Cell0D_coordinate(1,v), 2) +
                pow(mesh.Cell0D_coordinate(2,u) - mesh.Cell0D_coordinate(2,v), 2)
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
    unsigned int vertice_in_analisi = vertice_finale;
    while (vertice_in_analisi != vertice_iniziale) {
        if (pred.find(vertice_in_analisi) == pred.end()) {
            cout << "Nessun cammino trovato!" << endl;
            return;
        }
        path.push_back(vertice_in_analisi);
        vertice_in_analisi = pred[vertice_in_analisi];
    }
    path.push_back(vertice_iniziale);
    reverse(path.begin(), path.end());

    // Salva risultati nella mesh
    mesh.num_archiPath = path.size() - 1;
    mesh.lunghezza_Path = dist[vertice_finale];

    // Marca i vertici nel cammino
    mesh.Cell0D_marker[1].clear();
    for (unsigned int v : path) {
        mesh.Cell0D_marker[1].push_back(v);
        cout << "Vertice cammino: " << v << endl;
    }

    // Marca gli archi nel cammino
    mesh.Cell1D_marker[1].clear();
    for (size_t i = 1; i < path.size(); ++i) {
        unsigned int u = path[i-1];
        unsigned int v = path[i];
        for (unsigned int j = 0; j < mesh.Cell1D_num; ++j) {
            if ((mesh.Cell1D_estremi(0, j) == u && mesh.Cell1D_estremi(1, j) == v) ||
                (mesh.Cell1D_estremi(0, j) == v && mesh.Cell1D_estremi(1, j) == u)) {
                mesh.Cell1D_marker[1].push_back(mesh.Cell1D_id[j]);
                cout << "Arco cammino: " << j << endl;
                break;
            }
        }
    }
}
