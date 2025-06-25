#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <gtest/gtest.h>

#include "Utils.hpp"

/*-----------------------------------------------------------------------------------------------*/

// Test di integrità della struttura base PolyhedralMesh
TEST(PolyhedralMeshTest, StructureIntegrity) {
    PolyhedralLibrary::PolyhedralMesh mesh;
    mesh.Cell0D_num = 2;
    mesh.Cell0D_id = {0, 1};
    mesh.Cell0D_coordinate.resize(3, 2);
    mesh.Cell0D_coordinate.col(0) << 0.0, 0.0, 0.0;
    mesh.Cell0D_coordinate.col(1) << 1.0, 0.0, 0.0;

    mesh.Cell1D_num = 1;
    mesh.Cell1D_id = {0};
    mesh.Cell1D_estremi.resize(2, 1);
    mesh.Cell1D_estremi.col(0) << 0, 1;

    EXPECT_EQ(mesh.Cell0D_coordinate.cols(), mesh.Cell0D_num);
    EXPECT_EQ(mesh.Cell1D_estremi.cols(), mesh.Cell1D_num);
}

// Test costruzione ottaedro con triangolazione semplice
TEST(BuildSolidoTest, OctahedronTriangulatedOnce) {
    PolyhedralLibrary::PolyhedralMesh mesh;
    EXPECT_NO_THROW(build_classe_1(3, 4, 1, 0, mesh));

    EXPECT_EQ(mesh.Cell0D_num, 6);
    EXPECT_EQ(mesh.Cell1D_num, 12);
    EXPECT_EQ(mesh.Cell2D_num, 8);

    for (unsigned int i = 0; i < mesh.Cell0D_num; ++i)
        EXPECT_EQ(mesh.Cell0D_id[i], i);
    for (unsigned int i = 0; i < mesh.Cell1D_num; ++i)
        EXPECT_EQ(mesh.Cell1D_id[i], i);
    for (unsigned int i = 0; i < mesh.Cell2D_num; ++i)
        EXPECT_EQ(mesh.Cell2D_id[i], i);
}

// Test coerenza salvataggio file
TEST(BuildSolidoTest, FilesConsistencyWithMesh) {
    PolyhedralLibrary::PolyhedralMesh mesh;
    build_classe_1(3, 4, 1, 0, mesh);

    genera_file_Cell0D(mesh, "Cell0D.txt");
    genera_file_Cell1D(mesh, "Cell1D.txt");
    genera_file_Cell2D(mesh, "Cell2D.txt");

    std::ifstream file0D("Cell0D.txt");
    std::ifstream file1D("Cell1D.txt");
    std::ifstream file2D("Cell2D.txt");

    int count0D = std::count(std::istreambuf_iterator<char>(file0D), {}, '\n') - 1;
    int count1D = std::count(std::istreambuf_iterator<char>(file1D), {}, '\n') - 1;
    int count2D = std::count(std::istreambuf_iterator<char>(file2D), {}, '\n') - 1;

    EXPECT_EQ(count0D, mesh.Cell0D_num);
    EXPECT_EQ(count1D, mesh.Cell1D_num);
    EXPECT_EQ(count2D, mesh.Cell2D_num);
}

// Test Dijkstra su triangolo
TEST(DijkstraTest, ShortestPath_UsesDefinedFaceEdges) {
    using namespace PolyhedralLibrary;
    PolyhedralMesh mesh;

    mesh.Cell0D_num = 3;
    mesh.Cell0D_id = {0, 1, 2};
    mesh.Cell0D_coordinate.resize(3, 3);
    mesh.Cell0D_coordinate << 
        0.0, 1.0, 0.5,                     // vertice 0
        0.0, 0.0, std::sqrt(3.0)/2.0,                     // vertice 1
        0.0, 0.0, 0.0;      // vertice 2

    mesh.Cell1D_num = 3;
    mesh.Cell1D_id = {0, 1, 2};
    mesh.Cell1D_estremi.resize(2, 3);
    mesh.Cell1D_estremi << 0, 1, 2, 1, 2, 0;

    mesh.Cell2D_num = 1;
    mesh.Cell2D_vertici = {{0, 1, 2}};

    Dijkstra(mesh, 0, 1);

    EXPECT_EQ(mesh.num_archiPath, 1);
    EXPECT_NEAR(mesh.lunghezza_Path, 1.00000, 1e-15);
    ASSERT_TRUE(mesh.Cell0D_marker.count(1));
    EXPECT_EQ(mesh.Cell0D_marker[1].size(), 2);
    ASSERT_TRUE(mesh.Cell1D_marker.count(1));
    EXPECT_EQ(mesh.Cell1D_marker[1].size(), 1);
}

// Test duale: ottaedro -> cubo
TEST(DualeTest, OctahedronToCube) {
    using namespace PolyhedralLibrary;
    PolyhedralMesh mesh;

    std::vector<std::array<double, 3>> vertici = {
        {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}
    };
    mesh.Cell0D_num = vertici.size();
    mesh.Cell0D_id.resize(mesh.Cell0D_num);
    mesh.Cell0D_coordinate.resize(3, mesh.Cell0D_num);
    for (unsigned int i = 0; i < mesh.Cell0D_num; ++i) {
        mesh.Cell0D_id[i] = i;
        mesh.Cell0D_coordinate(0, i) = vertici[i][0];
        mesh.Cell0D_coordinate(1, i) = vertici[i][1];
        mesh.Cell0D_coordinate(2, i) = vertici[i][2];
    }

    std::vector<std::array<int, 3>> facce = {
        {0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
        {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}
    };
    mesh.Cell2D_num = facce.size();
    mesh.Cell2D_id.resize(mesh.Cell2D_num);
    mesh.Cell2D_vertici.resize(mesh.Cell2D_num);
    for (unsigned int i = 0; i < mesh.Cell2D_num; ++i) {
        mesh.Cell2D_id[i] = i;
        mesh.Cell2D_vertici[i] = {facce[i][0], facce[i][1], facce[i][2]};
    }

    PolyhedralMesh duale;
    build_duale(mesh, duale);

    EXPECT_EQ(duale.Cell0D_num, 8);
    EXPECT_EQ(duale.Cell1D_num, 12);
    EXPECT_EQ(duale.Cell2D_num, 6);

    for (unsigned int i = 0; i < duale.Cell0D_num; ++i) {
        double x = duale.Cell0D_coordinate(0, i);
        double y = duale.Cell0D_coordinate(1, i);
        double z = duale.Cell0D_coordinate(2, i);
        double norm = std::sqrt(x*x + y*y + z*z);
        EXPECT_NEAR(norm, 1.0, 1e-6);
    }
}

// Test classe 2 con b=1
TEST(BuildClasse2Test, OctahedronB1_PreciseCounts) {
    using namespace PolyhedralLibrary;
    PolyhedralMesh mesh;
    build_classe_2(3, 4, 1, 1, mesh); 

    EXPECT_EQ(mesh.Cell0D_num, 14);
    EXPECT_EQ(mesh.Cell1D_num, 36);
    EXPECT_EQ(mesh.Cell2D_num, 8);

    for (unsigned int i = 0; i < mesh.Cell0D_num; ++i)
        EXPECT_EQ(mesh.Cell0D_id[i], i);
    for (unsigned int i = 0; i < mesh.Cell1D_num; ++i)
        EXPECT_EQ(mesh.Cell1D_id[i], i);
    for (unsigned int i = 0; i < mesh.Cell2D_num; ++i)
        EXPECT_EQ(mesh.Cell2D_id[i], i);
}
