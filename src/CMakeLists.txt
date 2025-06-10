list(APPEND polyhedralMesh_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)
list(APPEND polyhedralMesh_sources ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp)

list(APPEND polyhedralMesh_headers ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedralMesh.hpp)

list(APPEND polyhedralMesh_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(polyhedralMesh_sources ${polyhedralMesh_sources} PARENT_SCOPE)
set(polyhedralMesh_headers ${polyhedralMesh_headers} PARENT_SCOPE)
set(polyhedralMesh_includes ${polyhedralMesh_includes} PARENT_SCOPE)
