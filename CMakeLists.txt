cmake_minimum_required(VERSION 2.8.3)
project(SparseReconstruction)

#add_custom_target(
#  DEPS ALL
#  COMMAND echo "downloading dependencies"
#  COMMAND cd ${PROJECT_SOURCE_DIR}/maxflow && pwd && make download_source
#)

ADD_SUBDIRECTORY (maxflow)

ADD_LIBRARY( SPARSE_RECONSTRUCTION STATIC
       src/SparseReconstruction.cpp
       include/SparseReconstruction.h
       )

include_directories(
  include
  #${PROJECT_SOURCE_DIR}/maxflow
) 

SET(GCC_CGAL_COMPILE_FLAGS "-frounding-math")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_CGAL_COMPILE_FLAGS}" )


set(SPARSE_RECONSTRUCTION_LIBRARIES GL glut gmp CGAL CGAL_Core)

TARGET_LINK_LIBRARIES( SPARSE_RECONSTRUCTION LINK_PUBLIC maxflow SPARSE_RECONSTRUCTION_LIBRARIES)
