
//-*- SparseReconstruction Class -*-
// Copyright 2016 Seyed Abbas Sadat (asadat@autonomylab.org)

#ifndef __SPARSE_RECONSTRUCTION_H
#define __SPARSE_RECONSTRUCTION_H

#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include "graph.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/metric_tsp_approx.hpp>

class Map;
class MapPoint;

using namespace TooN;

class SparseReconstruction
{
private:

    /*
     * stores a triangle and properties that indicate confidence
     */
    struct Triangle3D;

    public:

        class MapPoint_KeyFrames
        {
        public:
            // the 3D position of the map point
            Vector<3> p;

            // the 3D position of the viewpoints that
            // the map point was observed from
            std::vector<Vector<3> > keyframes_p;
        };

        class CellInfo{
        public:
            int id;
            static const int NO_ID = -1;
            std::set<int> collisionData;// the neighbor cells with which
                                           // this cell makes a surface boundary
            std::vector< Triangle3D* > surfaceTris; //the surface triangles that are on a face of this cell.
            CellInfo();

            //returns the surface triangle shared by this cell and the cell with the supplied ID.
            //A NULL pointer is returned if the two cells don't share a surface triangle.
            Triangle3D* getCommonTriangle3D( int cell_id );
        };


        unsigned int nextCellId; //a global counter for assigning cell IDs

        class VertexInfo{
        private:
            static unsigned long nextVertexId; //a global counter for assigning vertex IDs
        public:
            long id;
            static const int NO_ID = -1;
            std::vector< Vector<3> > camera_locations; //the camera positions that saw this vertex
            std::vector< Triangle3D* > surfaceTris; //the surface triangles that have this vertex as a vertex
            VertexInfo();
        };


        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef CGAL::Triangulation_cell_base_with_info_3<CellInfo,K>       InfoCell;
        typedef CGAL::Triangulation_vertex_base_with_info_3<VertexInfo, K>  InfoVert;

        typedef CGAL::Triangulation_data_structure_3<InfoVert, InfoCell>         Tds;
        typedef CGAL::Delaunay_triangulation_3<K, Tds>                      Delaunay;

        typedef Delaunay::Point   Point;

        //this is the graph type used by the maxflow library
        typedef Graph<double,double,double> GraphType;


        ~SparseReconstruction();
        static SparseReconstruction * GetInstance()
        {
            if(instance == NULL)
            {
                instance = new SparseReconstruction();
            }

            return instance;
        }


        void glDraw();
        void Triangulate(std::vector<MapPoint_KeyFrames*> &mp_kf);

        void Reset();

        //bool saveMesh(const std::string& fileName);

        //bool savePoints(const std::string& file);

        //set the maximum distance a triangle can be from the collision triangle and still be considered a part of the surface estimation
        static const double SURFACE_EST_MAX_DIST = 1.0;

    private:

        struct Triangle3D{

            int id;
            Delaunay::Cell_handle c1,c2;
            int KeyframeTex;
            Vector<2> textureCord[3];
            Matrix<3,3> tri;
            std::vector<Delaunay::Vertex_handle> vertices;
            Delaunay::Vertex_handle outside_vh; //vertex not on this triangle, but on the cell in the OUTSIDE direction
            bool too_big;
            bool too_far;
            bool occluded; //assessed by IsSurfaceTriangleValid
            bool unexplored;
            Vector<3> unit_normal_towards_outside;

            //surface estimation fields
            Vector<3> estimate_plane_norm;
            Vector<3> estimate_plane_center;
            double estimate_plane_density;
            std::vector< Matrix<3,3> > neighbor_tris;
        };

        SparseReconstruction();
        static SparseReconstruction* instance;

        void GenerateGraph();

        void DrawCell(Delaunay::Cell_handle ch);
        Delaunay::Cell_handle GetCameraCell(const Vector<3> start, const Vector<3> end, Vector<3> & cpoint);

        Delaunay::Cell_handle FindCollisionFromVertex(const Vector<3> start, const Vector<3> end, const Delaunay::Cell_handle & cell,
                                                      Vector<3> & collisionPoint, Delaunay::Cell_handle & lastCell);

        Delaunay::Cell_handle NextCell(const Vector<3> start, const Vector<3> end, const Delaunay::Cell_handle & curCell,
                                       Vector<3> & cpoint1, Vector<3> & cpoint2);

        bool ExistsVertex(const Delaunay::Cell_handle &ch, const Vector<3> v);
        void copy(Vector<3> &v, SparseReconstruction::Delaunay::Point& p);
        int AssignIdToCell(const Delaunay::Cell_handle& cell);
        void IncreaseCost(const Delaunay::Cell_handle &cell1, const Delaunay::Cell_handle &cell2, double resetVal=-1);
        void IncreaseCost(int cell1, int cell2, double resetVal=-1);
        void AddSourceLink(int cellId);
        void AddSinkLink(int cellId);
        int GetGraphNodesCount() const;
        bool GetCommonTriangle(const Delaunay::Cell_handle &cell1, const Delaunay::Cell_handle &cell2,
                               Matrix<3,3> &triangle,bool noCameraNeighbor,
                               std::vector<Delaunay::Vertex_handle> &vertices, Delaunay::Vertex_handle &outside_vh);

        double MinDistToVertex(Delaunay & t, Point p);

        TooN::Vector<3> sf;

        Delaunay delaunayMesh;

        Vector<3> rstart;
        Vector<3> rend;
        std::vector<Delaunay::Cell_handle> collidedCells;

        double maxCost;
        std::map<std::pair<int,int>, double> edgeCost;
        std::map<int,Delaunay::Cell_handle> id2cell;
        std::map<int, std::pair<double, double> > sourceSinkCost;
        std::vector<Delaunay::Vertex_handle> cameraVertices;
        //Surface Triangle
        std::vector< Triangle3D* > surfaceTriangles;

        //min cut result
        std::vector< std::pair<Delaunay::Cell_handle,Delaunay::Cell_handle> > neighbor_cells;
        std::vector< Delaunay::Cell_handle > sourceNeighbors;
        std::vector< Delaunay::Cell_handle > sinkNeighbors;

        //graph functions
        void get_links_in_min_cut(int num_nodes, std::map< std::pair<int,int>, double > linksToCosts,
                                                   std::map<int, std::pair<double,double> > terminal_costs,
                                                   std::vector< std::pair<int,int> > &neighbor_links,
                                                   std::vector< int > &source_links_in_min_cut,
                                                   std::vector< int > &sink_links_in_min_cut);
        GraphType::node_id add_node_if_necessary(int n, std::map<int, GraphType::node_id> &node_id_map, GraphType *g);

        typedef boost::adjacency_matrix<boost::undirectedS, boost::no_property,
                boost::property <boost::edge_weight_t, double,
                boost::property<boost::edge_index_t, int> > > Boost_Graph;



        Vector<3> cgal2ToonVect( K::Vector_3 vect);

        inline Vector<3> cgalVect2ToonVect( K::Vector_3 vect){
            return makeVector( vect[0], vect[1], vect[2]);
        }

        inline Delaunay::Point toonVect2cgalPoint( Vector<3> vect){
            return Delaunay::Point( vect[0], vect[1], vect[2]);
        }

        inline Vector<3> cgalPoint2ToonVect( K::Point_3 p){
            return makeVector( p[0], p[1], p[2]);
        }
        inline K::Triangle_3 toonMat2CgalTri( Matrix<3,3> tri_mat){
            K::Point_3 p1 = K::Point_3(tri_mat[0][0], tri_mat[0][1], tri_mat[0][2]);
            K::Point_3 p2 = K::Point_3(tri_mat[1][0], tri_mat[1][1], tri_mat[1][2]);
            K::Point_3 p3 = K::Point_3(tri_mat[2][0], tri_mat[2][1], tri_mat[2][2]);
            return K::Triangle_3(p1,p2,p3);
        }
        inline Vector<3> get_triangle_centroid( Matrix<3,3> tri_mat ){
            return makeVector( (tri_mat[0][0]+tri_mat[1][0]+tri_mat[2][0])/3,
                               (tri_mat[0][1]+tri_mat[1][1]+tri_mat[2][1])/3,
                               (tri_mat[0][2]+tri_mat[1][2]+tri_mat[2][2])/3);
        }
        inline double euclidean_distance_squared(Vector<3> first, Vector<3> second){
            Vector<3> diff = first - second;
            return diff*diff;
        }

        //for getting shortest path through triangles
        template<typename VertexListGraph, typename Vertex>
        std::map<Vertex, Triangle3D> map_vertices_to_triangles(std::vector<Triangle3D> tris, VertexListGraph &g);

        //adds edges between all vetrices having weights equal to the eclidean distance between the centers
        //of the triangles represented by the vertices.
        template<typename WeightMap, typename Vertex, typename VertexListGraph>
        void create_connected_graph(VertexListGraph &g, WeightMap wmap, std::map<Vertex,Triangle3D> vmap);


        Vector<3> get_norm_toward_outside(Triangle3D* tri);

        int intersect3D_RayTriangle( Vector<3> r0, Vector<3> r1 , Vector<3> v0, Vector<3> v1, Vector<3> v2, Vector<3> &I);
};



#endif
