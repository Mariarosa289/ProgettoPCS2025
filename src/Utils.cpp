#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace PolygonalLibrary
{
bool ImportMesh(PolygonalMesh& mesh)
{

    if(!ImportCell0Ds(mesh))
        return false;

    if(!ImportCell1Ds(mesh))
        return false;

    if(!ImportCell2Ds(mesh))
        return false;

    return true;

}
// ***************************************************************************
bool ImportCell0Ds(PolygonalMesh& mesh)
{
	

    ifstream file("./Cell0Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;

    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell0Ds = listLines.size();

    if (mesh.NumCell0Ds == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }

    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);
		string token;
		

        unsigned int id;
        unsigned int marker;
        Vector2d coord;
		
		getline(converter, token, ';');
		id = stoi(token);

		getline(converter, token, ';');
		marker = stoi(token);

		getline(converter, token, ';');
		mesh.Cell0DsCoordinates(0, id) = stod(token);
		
		getline(converter, token, ';');
		mesh.Cell0DsCoordinates(1, id) = stod(token);

       

        mesh.Cell0DsId.push_back(id);

        /// Memorizza i marker
        if(marker != 0)
        {
            const auto it = mesh.MarkerCell0Ds.find(marker);
            if(it == mesh.MarkerCell0Ds.end())
            {
                mesh.MarkerCell0Ds.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell0Ds[marker].push_back(id);
                it->second.push_back(id);
            }
        }

    }

    return true;
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
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds);

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
