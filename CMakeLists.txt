set(CMAKE_BUILD_TYPE Release)

cmake_minimum_required(VERSION 3.9)

project(PolyhedralMesh LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
################################################################################

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest
find_package(GTest REQUIRED)

list(APPEND PolyhedralMesh_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND PolyhedralMesh_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

# Insert Sources
################################################################################
add_subdirectory(src)

list(APPEND PolyhedralMesh_SOURCES ${polyhedralMesh_sources})
list(APPEND PolyhedralMesh_HEADERS ${polyhedralMesh_headers})
list(APPEND PolyhedralMesh_INCLUDE ${polyhedralMesh_includes})

add_subdirectory(ExportParaview)

list(APPEND PolyhedralMesh_SOURCES ${src_paraview_sources})
list(APPEND PolyhedralMesh_HEADERS ${src_paraview_headers})
list(APPEND PolyhedralMesh_INCLUDE ${src_paraview_includes})


# Create executable
################################################################################
add_executable(${PROJECT_NAME}
      main.cpp
      ${PolyhedralMesh_SOURCES}
      ${PolyhedralMesh_HEADERS})


target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

add_executable(${PROJECT_NAME}_TEST
    main_test.cpp
    ${PolyhedralMesh_SOURCES}
    ${PolyhedralMesh_HEADERS}
    ${PolyhedrallMesh_TEST_SOURCES}
    ${PolyhedralMesh_HEADERS})

target_link_libraries(${PROJECT_NAME}_TEST ${PolyhedralMesh_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME}_TEST PRIVATE src_test ${PolyhedralMesh_INCLUDE} ${PolyhedralMesh_TEST_INCLUDE})
target_compile_options(${PROJECT_NAME}_TEST PUBLIC -fPIC)

target_link_libraries(${PROJECT_NAME} ${PolyhedralMesh_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${PolyhedralMesh_INCLUDE})
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)



