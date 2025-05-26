//RIVEDERE le #include

#include "Utils.hpp"
#include "PolygonalMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <tuple>

using namespace std;

//***************************************************************************
/// Funzione per normalizzare un punto su sfera = OK

array<double, 3> Normalizza(const array<double, 3>& p) {
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    return {p[0]/norm, p[1]/norm, p[2]/norm};
}

//***************************************************************************
/// Controllo input: false se input non validi = OK

bool ControllaInput(unsigned int p, unsigned int q, unsigned int b, unsigned int c) {
    return p >= 3 && q >= 3 && q <= 5 && (b == 0 || c == 0 || b == c);
}

//***************************************************************************
/// Genera l'icosaedro iniziale NON NORMALIZZATO centrato nell'origine

void build_ico(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce) {
    const double phi = (1.0 + sqrt(5.0)) / 2.0;   // per calcolare la sezione aurea phi => per avere l'icosaedro centrato in (0,0,0)
    vertici = {   //tutte combinazioni
        {-1,  phi,  0}, { 1,  phi,  0}, {-1, -phi,  0}, { 1, -phi,  0},
        { 0, -1,  phi}, { 0,  1,  phi}, { 0, -1, -phi}, { 0,  1, -phi},
        { phi,  0, -1}, { phi,  0,  1}, {-phi,  0, -1}, {-phi,  0,  1}
    };
    facce = {   //facce indicate con l'id dei vertici
        {0,11,5}, {0,5,1}, {0,1,7}, {0,7,10}, {0,10,11},
        {1,5,9}, {5,11,4}, {11,10,2}, {10,7,6}, {7,1,8},
        {3,9,4}, {3,4,2}, {3,2,6}, {3,6,8}, {3,8,9},
        {4,9,5}, {2,4,11}, {6,2,10}, {8,6,7}, {9,8,1}
    };
}

/// Genera l'ottaedro iniziale NORMALIZZATO centrato nell'origine (per proprietà dell'ottaedro)

void build_octa(vector<array<double, 3>>& vertici,vector<array<int, 3>>& facce) {
    vertici = {
        { 1,  0,  0},  // vertice 0
        {-1,  0,  0},  // vertice 1
        { 0,  1,  0},  // vertice 2
        { 0, -1,  0},  // vertice 3
        { 0,  0,  1},  // vertice 4 (alto)
        { 0,  0, -1}   // vertice 5 (basso)
    };

    // 8 facce triangolari, definite dai vertici
    facce = {
        {0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4}, // facce superiori
        {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}  // facce inferiori
    };
}

/// Genera il tetraedo iniziale NON NORMALIZZATO centrato nell'origine

void build_tetra(vector<array<double, 3>>& vertici, vector<array<int, 3>>& facce) {
    vertici = {
        {  1,  1,  1 },
        { -1, -1,  1 },
        { -1,  1, -1 },
        {  1, -1, -1 }
    };

    facce = {
        {0, 1, 2},
        {0, 3, 1},
        {0, 2, 3},
        {1, 3, 2}
    };
}


//****************************************************
/// Funzione per salvare i vertici normalizzati = OK

unsigned int salva_vertice_norm(const array<double,3>& coord,
                           map<array<double,3>, unsigned int>& map0D,
                           PolygonalMesh& mesh,
                           unsigned int& vid) 
    {
        auto norm = Normalizza(coord);
        if (map0D.count(norm)) 
            return map0D[norm];  // se c'è il punto, non aggiorna il vid
        Vertex v{vid, norm};
        mesh.Cell0D.push_back(v);
        return map0D[norm] = vid++;  // vid utilizzato -> aggiorniamo il vid
    }


// ***************************************************************************
/// Costruzione del geodetico classe I

void build_solido(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolygonalMesh& mesh) {   //prende gli input e costruisce il solido
    if (!ControllaInput(p, q, b, c)) throw runtime_error("Input non valido");   // SCRIVERE NEL MAIN
    if (p != 3 || c != 0) throw runtime_error("Supportati solo solidi di classe I (p=3, c=0)");   //SCRIVERE NEL MAIN

    vector<array<double, 3>> vertici_iniziali;
    vector<array<int, 3>> facce_iniziali;

    // costruisco vertici e facce iniziali
    if (q == 5) {
        build_ico(vertici_iniziali, facce_iniziali);
    } else if (q == 4) {
        build_octa(vertici_iniziali, facce_iniziali)  ; 
    } else {
        build_tetra(vertici_iniziali, facce_iniziali);   
    };


    mesh = PolygonalMesh();   //per riempire le CellxD


    unsigned int vid = 0, eid = 0, fid = 0;   //inizializzare gli ID dei vertici, dei lati e delle facce

    /// Mappa per Cell0D
    map<array<double, 3>, unsigned int> map0D;   //Dizionario che lega id e vertice

    for (const auto& face : facce_iniziali) {
        const auto& A = vertici_iniziali[face[0]];
        const auto& B = vertici_iniziali[face[1]];
        const auto& C = vertici_iniziali[face[2]];

        vector<vector<unsigned int>> griglia(b+1);
        for (unsigned int i = 0; i <= b; ++i) {
            griglia[i].resize(i+1);
            for (unsigned int j = 0; j <= i; ++j) {
                double u = 1.0 - static_cast<double>(i)/b; // coordinate baricentriche 
                double v = static_cast<double>(i - j)/b; // static_cast <double> serve per far si che la divisione dia un risultato double
                double w = static_cast<double>(j)/b;
                array<double,3> P = {
                    u*A[0] + v*B[0] + w*C[0],
                    u*A[1] + v*B[1] + w*C[1],
                    u*A[2] + v*B[2] + w*C[2]
                };
                griglia[i][j] = salva_vertice_norm(P, map0D, mesh, vid);  
            }
        }
// TRIANGOLAZIONE DELLA FACCIA PRINCIPALE IN BASE A b (salvando ogni sottotriangolo/faccia attraverso i vertici)
        for (unsigned int i = 0; i < b; ++i) { //itera sui vertici  dei sottotriangoli attraverso lo schema sugli appunti
            for (unsigned int j = 0; j < i + 1; ++j) {
                unsigned int v1 = griglia[i][j];
                unsigned int v2 = griglia[i+1][j];
                unsigned int v3 = griglia[i+1][j+1];
                mesh.Cell2D.push_back({fid++, {v1,v2,v3}, {}});
                if (j < i) {
                    unsigned int v4 = griglia[i][j+1];
                    mesh.Cell2D.push_back({fid++, {v1,v3,v4}, {}});
                }
            }
        }
    }
    /// Mappa per Cell1D
    map<pair<unsigned int,unsigned int>, unsigned int> map1D;
    for (auto& f : mesh.Cell2D) {
        for (int i = 0; i < 3; ++i) {
            unsigned int a = f.Cell0D[i];// serve per collegare tutti i vertici del triangolo
            unsigned int b = f.Cell0D[(i+1)%3];
            if (a > b) swap(a, b);// fa si che non ci siano duplicati essendo i lati non orientati
            auto key = make_pair(a, b); // crea la coppia senza specificare il tipo di a e di b
            if (!map1D.count(key)) {   // qui vede se c'è la key e se non c'è la aggiunge
                mesh.Cell1D.push_back({eid, a, b});
                map1D[key] = eid++;
            }
            f.Cell1D.push_back(map1D[key]);
        }
    }
/// COSTRUZIONE DEL POLIEDRO
    Cell3D poly{0}; // Inizializzo il poliedro
    for (auto& v : mesh.Cell0D) poly.Cell0D.push_back(v.id);
    for (auto& e : mesh.Cell1D) poly.Cell1D.push_back(e.id);
    for (auto& f : mesh.Cell2D) poly.Cell2D.push_back(f.id);
    mesh.Cell3D.push_back(poly);

}





