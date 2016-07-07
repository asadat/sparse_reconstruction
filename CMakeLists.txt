cmake_minimum_required(VERSION 2.8.3)
project(sparse_reconstruction)

ADD_DEFINITIONS(-DGL_VISUALIZATION)

ADD_LIBRARY( sparse_reconstruction SHARED
       maxflow/maxflow 
       maxflow/block.h 
       maxflow/graph.h 
       maxflow/instances.inc 
       maxflow/graph.cpp 
       maxflow/maxflow.cpp
       src/SparseReconstruction.cpp
       include/SparseReconstruction.h
       )

include_directories(
  include
  maxflow
) 



SET(GCC_CGAL_COMPILE_FLAGS "-frounding-math")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_CGAL_COMPILE_FLAGS} -Wno-write-strings -o3" )


set(SPARSE_RECONSTRUCTION_LIBRARIES GL glut gmp CGAL CGAL_Core)

set_target_properties(sparse_reconstruction
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
)

TARGET_LINK_LIBRARIES(sparse_reconstruction ${SPARSE_RECONSTRUCTION_LIBRARIES})

