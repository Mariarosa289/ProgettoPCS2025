#pragma once

#include <iostream>
#include <fstream>

#include <vector>

#include <gtest/gtest.h>

# include "Utils.hpp"

/*-----------------------------------------------------------------------------------------------*/

/// TEST per verificare il riempimento corretto della struct
	// inserimento di un singolo vertice in Cell0D
TEST(Cell0DTest, SingleVertexPopulation) { 
//Il primo parametro è il gruppo di test, e il secondo è il nome del singolo test.
    PolyhedralLibrary::PolyhedralMesh mesh;
	// inizializzo il numero di celle a 1
    mesh.NumCell0Ds = 1; 
	// il primo id di Cell0DsId sarà 0
    mesh.Cell0DsId.push_back(0); 
    mesh.Cell0DsCoordinates.resize(3,1);
	// coordinate del vertice di id=0
    mesh.Cell0DsCoordinates.col(0) << 1.0, 2.0, 3.0; 
	// controllo che Il numero di vertici sia effettivamente 1
    EXPECT_EQ(mesh.NumCell0Ds, 1); 
	// controllo che L'ID del primo vertice sia 0
    EXPECT_EQ(mesh.Cell0DsId[0], 0); 
	// controllo che Le coordinate siano esattamente quelle che ho messo
    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(0,0), 1.0); 
    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(1,0), 2.0); 
    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(2,0), 3.0); 
}

	// inserimento di più vertici in Cell0D
TEST(Cell0DTest, MultipleVerticesPopulation) {
    PolyhedralLibrary::PolyhedralMesh mesh;
	
	// inizializzo il numero di celle a 3
    mesh.NumCell0Ds = 3;
    mesh.Cell0DsId = {0, 1, 2};
    mesh.Cell0DsCoordinates.resize(3, 3);

    mesh.Cell0DsCoordinates.col(0) << 1.0, 2.0, 3.0; // Vertex 0
    mesh.Cell0DsCoordinates.col(1) << 4.0, 5.0, 6.0; // Vertex 1
    mesh.Cell0DsCoordinates.col(2) << 7.0, 8.0, 9.0; // Vertex 2
	
	// controllo che Il numero di vertici sia effettivamente 3
    EXPECT_EQ(mesh.NumCell0Ds, 3);
    EXPECT_EQ(mesh.Cell0DsId.size(), 3);
    
    for (unsigned int i = 0; i < 3; ++i) {
		// controllo che L'ID del primo vertice sia 0, il secondo 1 e il terzo 2
        EXPECT_EQ(mesh.Cell0DsId[i], i);
    }
	
	// controllo che Le coordinate siano esattamente quelle che ho messo
    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(0,0), 1.0); // alla prima riga prima colonna ho 1.0 che è la coordinata x del vertice 0
    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(1,0), 2.0);
    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(2,0), 3.0);

    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(0,1), 4.0);
    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(1,1), 5.0);
    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(2,1), 6.0);

    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(0,2), 7.0);
    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(1,2), 8.0);
    EXPECT_DOUBLE_EQ(mesh.Cell0DsCoordinates(2,2), 9.0);
}

	// inserimento di uno singolo spigolo in Cell1D

TEST(Cell1DTest, SingleEdgePopulation) {
    PolyhedralLibrary::PolyhedralMesh mesh;
	
	// inizializzo il numero di spigoli a 1
    mesh.NumCell1Ds = 1;
    mesh.Cell1DsId.push_back(0);
    mesh.Cell1DsExtrema.resize(2, 1);  // 2 righe (da, a ->extrema) x 1 colonna (un solo spigolo)
    mesh.Cell1DsExtrema.col(0) << 0, 1; // da 0 a 1
	
	// controllo che Il numero di spigoli sia effettivamente 1
    EXPECT_EQ(mesh.NumCell1Ds, 1);
	// controllo che L'ID del primo spigolo sia 0
    EXPECT_EQ(mesh.Cell1DsId[0], 0);
	// controllo che gli extrema siano esattamente quelle che ho messo
    EXPECT_EQ(mesh.Cell1DsExtrema(0, 0), 0);  // da
    EXPECT_EQ(mesh.Cell1DsExtrema(1, 0), 1);  // a
}

	// inserimento di più singoli in Cell1D

TEST(Cell1DTest, MultipleEdgePopulation) {
    PolyhedralLibrary::PolyhedralMesh mesh;

	// inizializzo il numero di spigoli a 1
    mesh.NumCell1Ds = 3;
    mesh.Cell1DsId = {0, 1, 2};
    mesh.Cell1DsExtrema.resize(2, 3); // 2 righe (da, a ->extrema) x 3 colonna (3 spigoli)
    mesh.Cell1DsExtrema << 0, 1, 2,   // da
                            1, 2, 3;   // a

    EXPECT_EQ(mesh.NumCell1Ds, 3);
    EXPECT_EQ(mesh.Cell1DsId.size(), 3);

	//// controllo che L'ID del primo spigolo sia 0,del secondo 1 e del terzo 3
	// controllo che gli extrema siano esattamente quelle che ho messo, gli extrem sono inseriti in modo che l'id di partenza corrisponda alla posizione del contatore e quello di arrivo alla posizione successiva
    for (unsigned int i = 0; i < 3; ++i) {
        EXPECT_EQ(mesh.Cell1DsId[i], i);
        EXPECT_EQ(mesh.Cell1DsExtrema(0, i), i);
        EXPECT_EQ(mesh.Cell1DsExtrema(1, i), i + 1);
    }
}

	// inserimento di una singola faccia in Cell2D
	// verifico il caso base del solido geodetico, in cui la triangolazione crea solo facce triangolari
TEST(Cell2DTest, SingleFacePopulation) {
    PolyhedralLibrary::PolyhedralMesh mesh;

    mesh.NumCell2Ds = 1;
    mesh.Cell2DsId.push_back(0);
    mesh.Cell2DsVertices = {{0, 1, 2}};  // id dei vertici del Triangolo
    mesh.Cell2DsNumEdges.push_back(3);	// 3 spigoli per ogni faccia
    mesh.Cell2DsEdges = {{0, 1, 2}};     // Id degli spigoli

	// verifico che il numero di facce sia 1
    EXPECT_EQ(mesh.NumCell2Ds, 1);
	// verifico che l'id della faccia sia zero
    EXPECT_EQ(mesh.Cell2DsId[0], 0);
	// verifco che il numero di vertici della faccia sia 3
    EXPECT_EQ(mesh.Cell2DsVertices[0].size(), 3);
	// verifico che il numero di lati della faccia sia tre
    EXPECT_EQ(mesh.Cell2DsNumEdges[0], 3);
	// verifico che la dimezione del vettore che identifica gli spigoli della faccia sia tre
    EXPECT_EQ(mesh.Cell2DsEdges[0].size(), 3);
}

	// inserimento di una più facce in Cell2D
	// verifico il caso di facce triangolari, quadrate ed esagonali, che sono quelle che si presentano nella costruzione del solido di base e del suo duale ( per la prima triangoolazione)

TEST(Cell2DTest, MultipleFacesPopulation) {
    PolyhedralLibrary::PolyhedralMesh mesh;

    mesh.NumCell2Ds = 3;
    mesh.Cell2DsId = {0, 1, 2};

    // Triangolo
    mesh.Cell2DsVertices.push_back({0, 1, 2});
    mesh.Cell2DsEdges.push_back({0, 1, 2});
    mesh.Cell2DsNumEdges.push_back(3);

    // Quadrato
    mesh.Cell2DsVertices.push_back({3, 4, 5, 6});
    mesh.Cell2DsEdges.push_back({3, 4, 5, 6});
    mesh.Cell2DsNumEdges.push_back(4);

    // Esagono
    mesh.Cell2DsVertices.push_back({7, 8, 9, 10, 11, 12});
    mesh.Cell2DsEdges.push_back({7, 8, 9, 10, 11, 12});
    mesh.Cell2DsNumEdges.push_back(6);

    // Controllo che il numero di facce inserite dia 3
    EXPECT_EQ(mesh.NumCell2Ds, 3);

    for (size_t i = 0; i < mesh.NumCell2Ds; ++i) {
        // ID coerente
        EXPECT_EQ(mesh.Cell2DsId[i], i);

        // Dati coerenti: la dimensione del vettore che contiene gli id dei lati di ogni faccia deve essere uguale al numero di lati di ogni faccia
        EXPECT_EQ(mesh.Cell2DsEdges[i].size(), mesh.Cell2DsNumEdges[i]);

        // Anche i vertici devono essere almeno tanti quanto gli spigoli
        EXPECT_GE(mesh.Cell2DsVertices[i].size(), mesh.Cell2DsNumEdges[i]);
    }
	
	
}


// testo se la costruzione dell'ottaedro avviene in maniera regolare, rispettando il numero di vertici, spigoli e facce 

TEST(BuildSolidoTest, OctahedronTriangulatedOnce) {
    PolyhedralLibrary::PolyhedralMesh mesh;

    // Ottaedro triangolato con 1 suddivisione
    EXPECT_NO_THROW(build_classe_1(3, 4, 1, 0, mesh));
	
    // Verifica esatta dei numeri
    EXPECT_EQ(mesh.NumCell0Ds, 6); // numero atteso di vertici
    EXPECT_EQ(mesh.NumCell1Ds, 12); // numero atteso di spigoli
    EXPECT_EQ(mesh.NumCell2Ds, 8); // numero atteso di facce

    // ID coerenti : gli id di vertici lati e facce devono essere numerati in maniera crescente da 0 a N-1
    for (unsigned int i = 0; i < mesh.NumCell0Ds; ++i)
        EXPECT_EQ(mesh.Cell0DsId[i], i);
    for (unsigned int i = 0; i < mesh.NumCell1Ds; ++i)
        EXPECT_EQ(mesh.Cell1DsId[i], i);
    for (unsigned int i = 0; i < mesh.NumCell2Ds; ++i)
        EXPECT_EQ(mesh.Cell2DsId[i], i);
	
	// testo la condizione faces.edges[e].end == faces.edges[(e + 1)%E].originfaces.vertices[e] == faces.edges[e].origin
	for (unsigned int f = 0; f < mesh.NumCell2Ds; ++f) {
    int E = mesh.Cell2DsNumEdges[f];
    
    for (int e = 0; e < E; ++e) {
        int edge_id = mesh.Cell2DsEdges[f][e];
        int next_edge_id = mesh.Cell2DsEdges[f][(e + 1) % E];
        
        // origin e end dello spigolo corrente
        int origin = mesh.Cell1DsExtrema(0, edge_id);
        int end = mesh.Cell1DsExtrema(1, edge_id);
        
        // origin dello spigolo successivo
        int next_origin = mesh.Cell1DsExtrema(0, next_edge_id);

        // Verifica faces.edges[e].end == faces.edges[(e + 1) % E].origin
        EXPECT_EQ(end, next_origin);

        // Verifica faces.vertices[e] == faces.edges[e].origin
        int vertex_id = mesh.Cell2DsVertices[f][e];
        EXPECT_EQ(vertex_id, origin);
    }
}



}





	// test input invalidi
	
	
/*TEST(BuildSolidoTest, InvalidInputs) {
    PolyhedralLibrary::PolyhedralMesh mesh;

    // p non valido
    EXPECT_THROW(wrapper(4, 4, 2, 0, mesh), std::runtime_error);

    // c non valido
    EXPECT_THROW(build_classe_1(3, 4, 2, 1, mesh), std::runtime_error);

    // ControllaInput fallisce
    EXPECT_THROW(build_classe_1(3, 0, 1, 0, mesh), std::runtime_error); 
}*/

// test per verificare costruzione coerente del file: controllo che il numero di righe nei file relativi ai vertici,spigoli e facce sia uguale al numero di vertci spigoli e facce del poliedro generato //

TEST(BuildSolidoTest, FilesConsistencyWithMesh) {
    PolyhedralLibrary::PolyhedralMesh mesh;
    build_classe_1(3, 4, 1, 0, mesh);

    // Genera i file, passandogli nome file
    GeneraFileCell0D(mesh, "Cell0D.txt");
    GeneraFileCell1D(mesh, "Cell1D.txt");
    GeneraFileCell2D(mesh, "Cell2D.txt");

    // Apri i file e conta le righe (assumendo una riga per elemento)
    std::ifstream file0D("Cell0D.txt");
    std::ifstream file1D("Cell1D.txt");
    std::ifstream file2D("Cell2D.txt");

	// conto le righe del file: conta quante volte compare il carattere \n (newline) nell’intervallo che va dall’inizio alla fine del file.
	// std::istreambuf_iterator<char>(file0D) iteratore che inizia a leggere il file dal primo carattere (quindi dall’inizio del file).
	// std::istreambuf_iterator<char>() iteratore  che indica la fine del file
	// tolgo uno perchè non va contato l'header
    int count0D = std::count(std::istreambuf_iterator<char>(file0D), std::istreambuf_iterator<char>(), '\n')-1;
    int count1D = std::count(std::istreambuf_iterator<char>(file1D), std::istreambuf_iterator<char>(), '\n')-1;
    int count2D = std::count(std::istreambuf_iterator<char>(file2D), std::istreambuf_iterator<char>(), '\n')-1;

    EXPECT_EQ(count0D, mesh.NumCell0Ds);
    EXPECT_EQ(count1D, mesh.NumCell1Ds);
    EXPECT_EQ(count2D, mesh.NumCell2Ds);
}

// TEST PER VERIFICARE LA CORRETTEZZA DELLA FUNZIONE PER IL CALCOLO DEL CAMMINO MINIMO//

TEST(DijkstraTest, ShortestPath_UsesDefinedFaceEdges) {
    using namespace PolyhedralLibrary;

    PolyhedralMesh mesh;

    // Definizione di 3 vertici di un triangolo equilatero
    mesh.NumCell0Ds = 3;
    mesh.Cell0DsId = {0, 1, 2};
    mesh.Cell0DsCoordinates.resize(3, 3);
    mesh.Cell0DsCoordinates << 
        0.0, 0.0, 0.0,                     // vertice 0
        1.0, 0.0, 0.0,                     // vertice 1
        0.5, std::sqrt(3.0)/2.0, 0.0;      // vertice 2

    // Spigoli ( ma servono per marcare il cammino)
    mesh.NumCell1Ds = 3;
    mesh.Cell1DsId = {0, 1, 2};
    mesh.Cell1DsExtrema.resize(2, 3);
    mesh.Cell1DsExtrema << 
        0, 1, 2,
        1, 2, 0; // Archi: 0-1, 1-2, 2-0

    // Una faccia con vertici (0, 1, 2) => crea archi: 0–1, 1–2, 2–0
    mesh.NumCell2Ds = 1;
    mesh.Cell2DsVertices = {{0, 1, 2}};

    // Esegui Dijkstra da vertice 0 a vertice 1
    Dijkstra(mesh, 0, 1);

    // Verifica che il cammino minimo sia diretto: 0 → 1 (distanza = 1.06488, non uno perchè i vertici vengono proiettati sulla sfera)
    EXPECT_EQ(mesh.num_archiPath, 1);
    EXPECT_NEAR(mesh.lunghezza_Path, 1.064882433048626, 1e-15);

    // Verifica che i vertici siano stati marcati correttamente
    ASSERT_TRUE(mesh.MarkerCell0Ds.count(1));
    EXPECT_EQ(mesh.MarkerCell0Ds[1].size(), 2);  // Vertici: 0, 1

    // Verifica che uno spigolo sia stato marcato
    ASSERT_TRUE(mesh.MarkerCell1Ds.count(1));
    EXPECT_EQ(mesh.MarkerCell1Ds[1].size(), 1);
}

// Test per costruzione corretta del duale

/*
TEST(DualeTest, OctahedronToCube) {
    using namespace PolyhedralLibrary;

    PolyhedralMesh mesh;

    // Vertici dell'ottaedro
    std::vector<std::array<double, 3>> vertici = {
        { 1,  0,  0}, {-1,  0,  0}, { 0,  1,  0},
        { 0, -1,  0}, { 0,  0,  1}, { 0,  0, -1}
    };
    mesh.NumCell0Ds = vertici.size();
    mesh.Cell0DsId.resize(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates.resize(3, mesh.NumCell0Ds);
    for (unsigned int i = 0; i < mesh.NumCell0Ds; ++i) {
        mesh.Cell0DsId[i] = i;
        mesh.Cell0DsCoordinates(0, i) = vertici[i][0];
        mesh.Cell0DsCoordinates(1, i) = vertici[i][1];
        mesh.Cell0DsCoordinates(2, i) = vertici[i][2];
    }

    // Facce dell'ottaedro (triangoli)
    std::vector<std::array<int, 3>> facce = {
        {0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
        {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}
    };
    mesh.NumCell2Ds = facce.size();
    mesh.Cell2DsId.resize(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.resize(mesh.NumCell2Ds);
    for (unsigned int i = 0; i < mesh.NumCell2Ds; ++i) {
        mesh.Cell2DsId[i] = i;
        mesh.Cell2DsVertices[i] = {facce[i][0], facce[i][1], facce[i][2]};
    }

    // Costruzione del duale
    PolyhedralMesh duale;
    build_duale(mesh, duale);

    // Verifiche : il duale dell'ottaedro è un cubo, che ha 8 vertici e 12 spigoli
    EXPECT_EQ(duale.NumCell0Ds, 8);  
    EXPECT_EQ(duale.NumCell1Ds, 12); 
    EXPECT_EQ(duale.NumCell2Ds, 6);  

    // Verifica che i vertici siano normalizzati (proiettati sulla sfera)
    for (unsigned int i = 0; i < duale.NumCell0Ds; ++i) {
        double x = duale.Cell0DsCoordinates(0, i);
        double y = duale.Cell0DsCoordinates(1, i);
        double z = duale.Cell0DsCoordinates(2, i);
        double norm = std::sqrt(x*x + y*y + z*z);
        EXPECT_NEAR(norm, 1.0, 1e-6) << "Vertice " << i << " non è normalizzato.";
    }

    


    // Verifica che tutti gli ID siano continui da 0 a N-1
    for (unsigned int i = 0; i < duale.NumCell0Ds; ++i)
        EXPECT_EQ(duale.Cell0DsId[i], i);
    for (unsigned int i = 0; i < duale.NumCell2Ds; ++i)
        EXPECT_EQ(duale.Cell2DsId[i], i);
}
        */













