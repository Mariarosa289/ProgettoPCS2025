#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <tuple>

using namespace std;

//***************************************************************************
// Funzione per normalizzare un punto su sfera = OK

array<double, 3> Normalize(const array<double, 3>& p) {
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    return {p[0]/norm, p[1]/norm, p[2]/norm};
}
//***************************************************************************
/// Controllo input: false se input non validi = OK

bool CheckInput(unsigned int p, unsigned int q, unsigned int b, unsigned int c) {
    return p >= 3 && q >= 3 && q <= 5 && (b == 0 || c == 0 || b == c);
}

//***************************************************************************
/// Genera l'icosaedro iniziale NON NORMALIZZATO = OK

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

//******************* */
/// Funzione per salvare i punti normalizzati : RIVEDERE PER Vertex

unsigned int salva_vertice(const array<double,3>& coord,
                           map<array<double,3>, unsigned int>& Cell0D,
                           PolygonalMesh& mesh,
                           unsigned int& vid) 
    {
        auto norm = Normalize(coord);
        if (Cell0D.count(norm)) 
            return Cell0D[norm];
        Vertex v{vid, norm};
        mesh.vertices.push_back(v);
        return Cell0D[norm] = vid++;
    }
// ***************************************************************************
/// Costruzione del geodetico classe I

void build_solido(unsigned int p, unsigned int q, unsigned int b, unsigned int c, PolygonalMesh& mesh) {   //prende gli input e costruisce il solido
    if (!CheckInput(p, q, b, c)) throw runtime_error("Input non valido");   // SCRIVERE NEL MAIN
    if (p != 3 || c != 0) throw runtime_error("Supportati solo solidi di classe I (p=3, c=0)");   //SCRIVERE NEL MAIN

    vector<array<double, 3>> vertici_iniziali;
    vector<array<int, 3>> facce_iniziali;

    // costruisco vertici e facce iniziali
    if (q == 5) {
        build_ico(vertici_iniziali, facce_iniziali)
    } else if (q == 4) {
        build_octa(vertici_iniziali, facce_iniziali)   //SCRIVERE CODICE DI build_octa
    } else {
        build_tetra(vertici_iniziali, facce_iniziali)   //SCRIVERE CODICE DI build_tetra
    };

    mesh = PolygonalMesh();


    unsigned int vid = 0, eid = 0, fid = 0;   //RIVEDERE

    // Inizializziamo Cell0D
    map<array<double, 3>, unsigned int> Cell0D;   //Dizionario che lega id e vertice



    for (const auto& face : facce_iniziali) {
        const auto& A = vertici_iniziali[face[0]];
        const auto& B = vertici_iniziali[face[1]];
        const auto& C = vertici_iniziali[face[2]];

        vector<vector<unsigned int>> grid(b+1);
        for (unsigned int i = 0; i <= b; ++i) {
            grid[i].resize(i+1);
            for (unsigned int j = 0; j <= i; ++j) {
                double u = 1.0 - static_cast<double>(i)/b;
                double v = static_cast<double>(i - j)/b;
                double w = static_cast<double>(j)/b;
                array<double,3> P = {
                    u*A[0] + v*B[0] + w*C[0],
                    u*A[1] + v*B[1] + w*C[1],
                    u*A[2] + v*B[2] + w*C[2]
                };
                grid[i][j] = salva_vertice(P);   //RIPARTIRE DA QUI
            }
        }

        for (unsigned int i = 0; i < b; ++i) {
            for (unsigned int j = 0; j < i + 1; ++j) {
                unsigned int v1 = grid[i][j];
                unsigned int v2 = grid[i+1][j];
                unsigned int v3 = grid[i+1][j+1];
                mesh.faces.push_back({fid++, {v1,v2,v3}, {}});
                if (j < i) {
                    unsigned int v4 = grid[i][j+1];
                    mesh.faces.push_back({fid++, {v1,v3,v4}, {}});
                }
            }
        }
    }

    map<pair<unsigned int,unsigned int>, unsigned int> edgeMap;
    for (auto& f : mesh.faces) {
        for (int i = 0; i < 3; ++i) {
            unsigned int a = f.vertices[i];
            unsigned int b = f.vertices[(i+1)%3];
            if (a > b) swap(a, b);
            auto key = make_pair(a, b);
            if (!edgeMap.count(key)) {
                mesh.edges.push_back({eid, a, b});
                edgeMap[key] = eid++;
            }
            f.edges.push_back(edgeMap[key]);
        }
    }

    Polyhedron poly{0};
    for (auto& v : mesh.vertices) poly.vertices.push_back(v.id);
    for (auto& e : mesh.edges) poly.edges.push_back(e.id);
    for (auto& f : mesh.faces) poly.faces.push_back(f.id);
    mesh.polyhedra.push_back(poly);

    ProjectVerticesOnSphere(mesh);
}





// ***************************************************************************
bool ImportCell1Ds(PolygonalMesh& mesh)
{
	cout << "ImportCell1Ds..." << endl;

    ifstream file("./Cell1Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();

    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::Matrixxi(2, mesh.NumCell1Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);
		string token;

        unsigned int id;
        unsigned int marker;
        Vector2i vertices;
		
		getline(converter, token, ';');
		id = stoi(token);

		getline(converter, token, ';');
		marker = stoi(token);

		getline(converter, token, ';');
		mesh.Cell1DsExtrema(0, id) = stod(token);
		
		getline(converter, token, ';');
		mesh.Cell1DsExtrema(1, id) = stod(token);
        
		cout << "id: " << id << ", marker: " << marker << ", extrema: (" 
		<< mesh.Cell1DsExtrema(0, id) << ", " 
		<< mesh.Cell1DsExtrema(1, id) << ")" << endl;
        mesh.Cell1DsId.push_back(id);

        /// Memorizza i marker
        if(marker != 0)
        {
            const auto it = mesh.MarkerCell1Ds.find(marker);
            if(it == mesh.MarkerCell1Ds.end())
            {
                mesh.MarkerCell1Ds.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell1Ds[marker].push_back(id);
                it->second.push_back(id);
            }
        }
    }
	for (const auto& entry : mesh.MarkerCell1Ds) {
		cout << "Marker: " << entry.first << " - IDs: ";
		for (const auto& id : entry.second)
			cout << id << " ";
		cout << endl;
	}
	cout << "Cell1DsExtrema dimensioni: " 
     << mesh.Cell1DsExtrema.rows() << "x" 
     << mesh.Cell1DsExtrema.cols() << endl;

    return true;
}
// ***************************************************************************
bool ImportCell2Ds(PolygonalMesh& mesh)
{
	cout << "ImportCell2Ds..." << endl;

    ifstream file;
    file.open("./Cell2Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell2Ds = listLines.size();

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);
		string token;
		
        unsigned int id;
		unsigned int marker;
		unsigned int NumVert;
		unsigned int NumEdg;
        vector<int> vertices;
        vector<int> edges;

        getline(converter, token, ';');
		id = stoi(token);

		getline(converter, token, ';');
		marker = stoi(token);

		getline(converter, token, ';');
		NumVert = stoi(token);
		
		vertices.resize(NumVert);
		
        for (unsigned int i = 0; i < NumVert; ++i)
		{
			getline(converter, token, ';');
			vertices[i] = stoi(token);
		}
		
		getline(converter, token, ';');
		NumEdg = stoi(token);
		edges.resize(NumEdg);
		
		
        for (unsigned int i = 0; i < NumEdg; ++i)
		{
			getline(converter, token, ';');
			edges[i] = stoi(token);
			if (edges[i] == 0)
			{
				return false;
			}
		}

        mesh.Cell2DsId.push_back(id);
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsEdges.push_back(edges);
		mesh.Cell2DsNumEdges.push_back(NumEdg);
				
		vertices.clear();
		edges.clear();
		/// Memorizza i marker
        if(marker != 0)
        {
            const auto it = mesh.MarkerCell2Ds.find(marker);
            if(it == mesh.MarkerCell2Ds.end())
            {
                mesh.MarkerCell2Ds.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell2Ds[marker].push_back(id);
                it->second.push_back(id);
            }
        }
    }
		
	

    return true;
}

bool Area(PolygonalMesh& mesh)

{
	unsigned int END = mesh.Cell2DsId.size();
	unsigned int STARTER=0;
	for(unsigned int i =0;i<END;++i)
	{
		unsigned int SizeArea=0;
		unsigned int EDGES=mesh.Cell2DsNumEdges[i];
		for (unsigned int j=STARTER;j<STARTER+EDGES-2;++j)
		{
			SizeArea+=mesh.Cell0DsCoordinates(0,j)*mesh.Cell0DsCoordinates(1,j+1)-mesh.Cell0DsCoordinates(0,j+1)*mesh.Cell0DsCoordinates(1,j);
			
		}
		
		SizeArea+=mesh.Cell0DsCoordinates(0,STARTER+EDGES-1)*mesh.Cell0DsCoordinates(1,STARTER)-mesh.Cell0DsCoordinates(0,STARTER)*mesh.Cell0DsCoordinates(1,STARTER+EDGES-1);
		STARTER+=EDGES;
		if (SizeArea == 0)
			return false;

		
		
		
	}
	return true;
}

}
