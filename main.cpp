#include <iostream>
#include <stdexcept>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
// #include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;

bool ControllaInput(unsigned int p, unsigned int q, unsigned int b, unsigned int c) {
    return p >= 3 && q >= 3 && q <= 5 && (b == 0 || c == 0 || b == c);
}
int main() {
	if (!ControllaInput(p, q, b, c)) throw runtime_error("Input non valido");   // SCRIVERE NEL MAIN
    if (p != 3 || c != 0) throw runtime_error("Supportati solo solidi di classe I (p=3, c=0)");   //SCRIVERE NEL MAIN
    try {
        // Parametri di input (puoi anche leggerli da terminale o file)
        unsigned int p = 3;   // Numero lati per faccia (solo 3 supportato)
        unsigned int q = 5;   // Numero facce per vertice (5 = icosaedro)
        unsigned int b = 3;   // Parametro di suddivisione
        unsigned int c = 0;   // Solo classe I (c=0)

        PolygonalMesh mesh;

        // Costruzione della mesh
        build_solido(p, q, b, c, mesh);

        // Output riepilogo
        cout << "Geodetico costruito con successo!" << endl;
        cout << "Numero vertici (Cell0D): " << mesh.Cell0D.size() << endl;
        cout << "Numero lati (Cell1D):    " << mesh.Cell1D.size() << endl;
        cout << "Numero facce (Cell2D):   " << mesh.Cell2D.size() << endl;
        cout << "Numero poliedri:         " << mesh.polyhedra.size() << endl;
    }
    catch (const exception& e) {
        cerr << "Errore: " << e.what() << endl;
        return 1;
    }

    return 0;
}

