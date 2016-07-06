
#include "SparseReconstruction.h"

// for linux dirs
//#include <unistd.h>
//#include <sys/stat.h>
//#include <sys/types.h>
//#include <dirent.h>

//#include <fstream>

#include <CGAL/Triangle_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Point_2.h>
#include <CGAL/Line_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/ch_melkman.h>
#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Dimension.h>

#include "GL/glut.h"

#define SOURCE_COST 10000
#define SINK_COST   0.1
#define SMALLVALUE  0.00001
#define SMALL_NUM   0.00000001
#define VIS_SCORE_RAYS  1

SparseReconstruction * SparseReconstruction::instance = NULL;

SparseReconstruction::SparseReconstruction()
{

}


SparseReconstruction::~SparseReconstruction()
{    
}

void SparseReconstruction::Reset()
{
    delaunayMesh.clear();
    cameraVertices.clear();
    id2cell.clear();

    neighbor_cells.clear();
    sourceNeighbors.clear();
    sinkNeighbors.clear();

    //collisionData.clear();
    for(unsigned int i=0; i<surfaceTriangles.size(); i++)
        delete surfaceTriangles[i];

    surfaceTriangles.clear();

}

void SparseReconstruction::Triangulate(std::vector<MapPoint_KeyFrames*> &mp_kf)
{
    delaunayMesh.clear();
    cameraVertices.clear();
    nextCellId = 0;

//    for(unsigned int j=0; j< map->vpKeyFrames.size();j++)
//    {
//        Vector<3> v = map->vpKeyFrames[j]->se3CfromW.inverse().get_translation();
//        Point p(v[0],v[1],v[2]);
//        cameraVertices.push_back(t.insert(p));
//    }

    //copy the map to the delaunay mesh
        // This becomes an "unsafe" copy in the sense that some points may not be correct due to
        // a bundle adjustment moving the points during the copy.  This will cause the mesh to
        // be slightly wrong for this triangulation.  However, bundle adjustments dont move points
        // very far, and so the error in the triangulation should be minimal.  The mesh will still
        // be correct in the sense that no part will disappear, but may have graphical bugs from
        // some points being moved while reading and not others.
        // to make it safe,
    // for each point in the map


    for(size_t i=0; i<mp_kf.size(); i++)
    {
        // get the position
        Vector<3> v = mp_kf[i]->p;
        // make a delanay point
        Point p(v[0],v[1],v[2]);

        //ROS_INFO("here1");
        double distSq = MinDistToVertex(delaunayMesh, p);
        if(distSq < 0.001)
            continue;

        //ROS_INFO("here2");

        // add it to the mesh and get the handle
        Delaunay::Vertex_handle vh = delaunayMesh.insert(p);

        std::copy(mp_kf[i]->keyframes_p.begin(), mp_kf[i]->keyframes_p.end(), std::back_inserter(vh->info().camera_locations));

    }



    GenerateGraph();

    id2cell.clear();

    Delaunay::Cell_iterator itci;
    for(itci = delaunayMesh.all_cells_begin(); itci != delaunayMesh.all_cells_end(); ++itci)
    {
        if ( itci->info().id != CellInfo::NO_ID ){
            id2cell[itci->info().id] = itci;
        }
    }

    neighbor_cells.clear();
    sourceNeighbors.clear();
    sinkNeighbors.clear();

    std::vector< std::pair<int,int> > neighbor_cells_i;
    std::vector< int > sourceNeighbors_i;
    std::vector< int > sinkNeighbors_i;

    this->get_links_in_min_cut(GetGraphNodesCount(),edgeCost,sourceSinkCost,neighbor_cells_i,sourceNeighbors_i,sinkNeighbors_i);

    for(unsigned int i=0; i<sinkNeighbors_i.size();i++)
    {
        sinkNeighbors.push_back(id2cell[sinkNeighbors_i[i]]);
    }

    for(unsigned int i=0; i<neighbor_cells_i.size(); i++)
    {
        std::pair<Delaunay::Cell_handle,Delaunay::Cell_handle> edge;

        edge.first  = id2cell[neighbor_cells_i[i].first];
        edge.second = id2cell[neighbor_cells_i[i].second];
        neighbor_cells.push_back(edge);

    }


    //collisionData.clear();
    for(unsigned int i=0; i<surfaceTriangles.size(); i++)
        delete surfaceTriangles[i];
    surfaceTriangles.clear();
    int tri_id = 0;
    for(unsigned int i=0; i < neighbor_cells.size(); i++)
    {
        Matrix<3,3> tr;
        Triangle3D *tri = new Triangle3D();

        if(GetCommonTriangle(neighbor_cells[i].first, neighbor_cells[i].second, tr, true, tri->vertices, tri->outside_vh))
        {
            //Add Collision data
            neighbor_cells[i].first->info().collisionData.insert(neighbor_cells[i].second->info().id);


            tri->id = tri_id;
            tri->tri = tr;
            tri->KeyframeTex = -1;//VisibilityServer::GetInstance()->GetBestKeyFrame(tr[0], tr[1], tr[2]);

            if(tri->KeyframeTex >= 0)
            {
                //tri->textureCord[0] = VisibilityServer::GetInstance()->GetImageCoord(tr[0], tri->KeyframeTex);
                //tri->textureCord[1] = VisibilityServer::GetInstance()->GetImageCoord(tr[1], tri->KeyframeTex);
                //tri->textureCord[2] = VisibilityServer::GetInstance()->GetImageCoord(tr[2], tri->KeyframeTex);
            }
            else
            {
                tri->textureCord[0] = makeVector(0,0);
                tri->textureCord[1] = makeVector(0,0);
                tri->textureCord[2] = makeVector(0,0);
            }

            tri->c1 = neighbor_cells[i].first;
            tri->c2 = neighbor_cells[i].second;

            tri->c1->info().surfaceTris.push_back(tri);
            tri->c2->info().surfaceTris.push_back(tri);

            //associate this triangle3d with each of its vertices
            for (std::vector<Delaunay::Vertex_handle>::iterator vit = tri->vertices.begin(); vit != tri->vertices.end(); vit++ )
            {
                (*vit)->info().surfaceTris.push_back( tri );
            }

            tri->unit_normal_towards_outside = get_norm_toward_outside( tri );
            TooN::normalize( tri->unit_normal_towards_outside );

            tri->unexplored = (tri->vertices[0]->info().camera_locations.size() <= 0) || (tri->vertices[1]->info().camera_locations.size() <= 0)
                    || (tri->vertices[2]->info().camera_locations.size() <= 0);

            assert( tri->outside_vh!=NULL ); //set by getCommonTriangle
            surfaceTriangles.push_back(tri);

            tri_id++;
        }
    }
}

void SparseReconstruction::copy(Vector<3> &v, Delaunay::Point &p)
{
    v[0] = p.x();
    v[1] = p.y();
    v[2] = p.z();
}

int SparseReconstruction::intersect3D_RayTriangle( Vector<3> r0, Vector<3> r1 , Vector<3> v0, Vector<3> v1, Vector<3> v2, Vector<3> &I)
{
    Vector<3>    u, v, n;              // triangle vectors
    Vector<3>    dir, w0, w;           // ray vectors
    float     r, a, b;              // params to calc ray-plane intersect

    // get triangle edge vectors and plane normal
    u = v1 - v0;
    v = v2 - v0;
    n = u ^ v;              // cross product
    if (n == (TooN::makeVector(0,0,0)))             // triangle is degenerate
        return -1;                  // do not deal with this case

    dir = r1 - r0;              // ray direction vector
    w0 = r0 - v0;
    a = -(n*w0);
    b = (n*dir);
    if (fabs(b) < SMALL_NUM) {     // ray is  parallel to triangle plane
        if (a == 0)                 // ray lies in triangle plane
            return 2;
        else return 0;              // ray disjoint from plane
    }

    // get intersect point of ray with triangle plane
    r = a / b;
    if (r < 0.0)                    // ray goes away from triangle
        return 0;                   // => no intersect
    // for a segment, also test if (r > 1.0) => no intersect

    I = r0 + r * dir;            // intersect point of ray and plane

    // is I inside T?
    float    uu, uv, vv, wu, wv, D;
    uu = (u*u);
    uv = (u*v);
    vv = (v*v);
    w = I - v0;
    wu = (w*u);
    wv = (w*v);
    D = (uv) *  uv - (uu) * vv;

    // get and test parametric coords
    float s, t;
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0)         // I is outside T
        return 0;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  // I is outside T
        return 0;

    return 1;                       // I is in T
}

Vector<3> SparseReconstruction::get_norm_toward_outside(Triangle3D *tri)
{
    //find normal in direction of OUTSIDE
    Vector<3> v0 = cgalPoint2ToonVect( tri->vertices[0]->point() );
    Vector<3> v1 = cgalPoint2ToonVect( tri->vertices[1]->point() );
    Vector<3> v2 = cgalPoint2ToonVect( tri->vertices[2]->point() );

    Vector<3> triNorm = (v0-v1)^(v0-v2);

    Vector<3> vOutside = cgalPoint2ToonVect( tri->outside_vh->point() );
    Vector<3> center =  0.333 * (v0 + v1 + v2);

    //choose triangle normal that points towards the outside vertex
    if(triNorm * (vOutside - center) < 0)
        triNorm = -1 * triNorm;

    return triNorm;
}

bool SparseReconstruction::ExistsVertex(const Delaunay::Cell_handle &ch, const Vector<3> v)
{
    Delaunay::Point p(v[0], v[1], v[2]);
    Delaunay::Vertex_handle vh;
    if(delaunayMesh.is_vertex(p, vh))
    {
        if(ch->has_vertex(vh))
            return true;
        else
            return false;
    }
    else
    {
        //ROS_INFO("Feature point not found in the mesh vertices.");
        return false;
    }
}

SparseReconstruction::Delaunay::Cell_handle SparseReconstruction::GetCameraCell(Vector<3> cameraPos, Vector<3> featurePose, Vector<3> & cpoint)
{
    Vector<3> near = featurePose-cameraPos;
    normalize(near);
    near = cameraPos + SMALLVALUE * near;
    cpoint = near;
    Delaunay::Point p(near[0],near[1],near[2]);
    Delaunay::Cell_handle cell = delaunayMesh.locate(p);
    return cell;
}

SparseReconstruction::Delaunay::Cell_handle SparseReconstruction::NextCell(const Vector<3> start, const Vector<3> end, const Delaunay::Cell_handle &curCell, Vector<3> &cpoint1, Vector<3> &cpoint2)
{
    Delaunay::Cell_handle ch;
    return ch;
}

SparseReconstruction::Delaunay::Cell_handle SparseReconstruction::FindCollisionFromVertex(const Vector<3> start, const Vector<3> end, const Delaunay::Cell_handle &cell,
                                                                                Vector<3> &cpoint, Delaunay::Cell_handle & lastCell)
{

    Delaunay::Vertex_handle vh[4];
    Vector<3> vp[4];
    for(int i=0; i<4; i++)
    {
        vh[i] = cell->vertex(i);
        copy(vp[i], vh[i]->point());
    }

    int res[] = {100,100,100,100};
    int checks[] = {0,0,0,0};
    //int ii=-1, jj=-1, kk=-1;

    Delaunay::Cell_handle nextCell = NULL;
    //conds[0] = (blackList.find(vh[0]) == blackList.end()) || (blackList.find(vh[1]) == blackList.end()) || (blackList.find(vh[2]) == blackList.end());
    ///conds[1] = (blackList.find(vh[0]) == blackList.end()) || (blackList.find(vh[1]) == blackList.end()) || (blackList.find(vh[3]) == blackList.end());
    ///conds[2] = (blackList.find(vh[0]) == blackList.end()) || (blackList.find(vh[2]) == blackList.end()) || (blackList.find(vh[3]) == blackList.end());
    ///conds[3] = (blackList.find(vh[1]) == blackList.end()) || (blackList.find(vh[2]) == blackList.end()) || (blackList.find(vh[3]) == blackList.end());

    bool foundCollision = false;



    {       
        res[0] = intersect3D_RayTriangle(start, end, vp[0], vp[1], vp[2], cpoint);

//        if(res[0] == 1)
//            collisionPoints.push_back(cpoint);

//        checks[0] += (cell->neighbor(3) != lastCell)?1:0;
//        checks[0] += (!ExistsVertex(cell,cpoint))?5:0;

        if(res[0] == 1 && cell->neighbor(3) != lastCell && !ExistsVertex(cell,cpoint))
        {
            nextCell = cell->neighbor(3);
//            ii=0; jj=1; kk=2;
//            blackList.clear();
//            blackList.insert(vh[0]);
//            blackList.insert(vh[1]);
//            blackList.insert(vh[2]);
            foundCollision = true;
        }
    }

    if(!foundCollision)
    {
        res[1] = intersect3D_RayTriangle(start, end, vp[0], vp[1], vp[3], cpoint);

//        if(res[1] == 1)
//            collisionPoints.push_back(cpoint);

//        checks[1] += (cell->neighbor(2) != lastCell)?1:0;
//        checks[1] += (!ExistsVertex(cell,cpoint))?5:0;

        if(res[1]==1 && cell->neighbor(2) != lastCell && !ExistsVertex(cell,cpoint))
        {
            nextCell = cell->neighbor(2);
//            ii=0; jj=1; kk=3;
//            blackList.clear();
//            blackList.insert(vh[0]);
//            blackList.insert(vh[1]);
//            blackList.insert(vh[3]);
            foundCollision = true;
        }

    }

    if(!foundCollision)
    {
        res[2] = intersect3D_RayTriangle(start, end, vp[0], vp[2], vp[3], cpoint);

//        if(res[2] == 1)
//            collisionPoints.push_back(cpoint);

//        checks[2] += (cell->neighbor(1) != lastCell)?1:0;
//        checks[2] += (!ExistsVertex(cell,cpoint))?5:0;

        if(res[2] == 1 && cell->neighbor(1) != lastCell && !ExistsVertex(cell,cpoint))
        {
            nextCell = cell->neighbor(1);
//            ii=0; jj=2; kk=3;
//            blackList.clear();
//            blackList.insert(vh[0]);
//            blackList.insert(vh[2]);
//            blackList.insert(vh[3]);
            foundCollision = true;
        }
    }

    if(!foundCollision)
    {
        res[3] = intersect3D_RayTriangle(start, end, vp[1], vp[2], vp[3], cpoint);

//        if(res[3] == 1)
//            collisionPoints.push_back(cpoint);

//        checks[3] += (cell->neighbor(0) != lastCell)?1:0;
//        checks[3] += (!ExistsVertex(cell,cpoint))?5:0;

        if(res[3] == 1 && cell->neighbor(0) != lastCell && !ExistsVertex(cell,cpoint))
        {
            nextCell = cell->neighbor(0);
//            ii=1; jj=2; kk=3;
//            blackList.clear();
//            blackList.insert(vh[1]);
//            blackList.insert(vh[2]);
//            blackList.insert(vh[3]);
            foundCollision = true;
        }
    }



/*
    if(ii>=0 && jj >=0 && kk>= 0)
    {
        for(int j=0; j<4; j++)
        {
            Delaunay::Cell_handle chn = cell->neighbor(j);
            if(chn->has_vertex(vh[ii]) && chn->has_vertex(vh[jj]) && chn->has_vertex(vh[kk]))
            {
                //if(!t.is_infinite(chn))
                {
                    nextCell = chn;

                }
//                else
//                {
//                    ROS_INFO("hitting infinite cell!!!!");
//                }
            }
        }
    }
    else
    {
        ROS_INFO("No Collision with the sides of he cell!!!! res: %d\t%d\t%d\t%d", res[0], res[1], res[2], res[3]);
    }
*/
    if(nextCell == NULL)
    {

        collidedCells.clear();
        collidedCells.push_back(lastCell);
        collidedCells.push_back(cell);
        rstart = start;
        rend = end;        
    }

    return nextCell;
}

int SparseReconstruction::AssignIdToCell(const Delaunay::Cell_handle &cell)
{
    if (cell->info().id == CellInfo::NO_ID)
    {
        cell->info().id  = nextCellId++;
    }
    return cell->info().id;

//    std::map<Delaunay::Cell_handle,int>::iterator it;
//    int s = cell2id.size();
//    //int curId
//    it = cell2id.find(cell);
//    if( it != cell2id.end())
//    {
//        return it->second;
//    }
//    else
//    {
//        cell2id[cell] = s;
//        return s;
//    }

  /*  bool flag = true;
    for(int i=0; i<collidedCells.size(); i++)
    {
        if(collidedCells[i] == cell)
        {
            flag = false;
            break;
        }
    }

    if(flag)
        collidedCells.push_back(cell);
        */
}

void SparseReconstruction::IncreaseCost(const Delaunay::Cell_handle &cell1, const Delaunay::Cell_handle &cell2, double resetVal)
{
    int cell1id = AssignIdToCell(cell1);
    int cell2id = AssignIdToCell(cell2);

    IncreaseCost(cell1id, cell2id, resetVal);
}

void SparseReconstruction::IncreaseCost(int cell1, int cell2, double resetVal)
{
    std::pair<int,int> edge;
    edge.first = cell1;
    edge.second = cell2;

    std::map<std::pair<int,int>, double>::iterator it;
    it = edgeCost.find(edge);
    double cost;

//    if(it == edgeCost.end())
//    {
//        edge.first = cell2;
//        edge.second = cell1;
//        it = edgeCost.find(edge);
//    }

    if(it != edgeCost.end())
    {
        if(resetVal>=0)
            it->second = resetVal;
        else
            it->second += 0.1;

        cost = it->second;
    }
    else
    {
        if(resetVal>=0)
            it->second = resetVal;
        else
            edgeCost[edge] = 0.1;

        cost = edgeCost[edge];
    }

    if(cost > maxCost && cost < SINK_COST)
        maxCost = cost;
}

void SparseReconstruction::AddSourceLink(int cellId)
{
    std::map<int, std::pair<double,double> >::iterator it;
    it = sourceSinkCost.find(cellId);
    if(it != sourceSinkCost.end())
    {
        it->second.first = SOURCE_COST;
    }
    else
    {
        std::pair<double,double> cost;
        cost.first = SOURCE_COST;
        cost.second = 0;
        //it->second = cost;
        sourceSinkCost[cellId] = cost;
    }

}

int SparseReconstruction::GetGraphNodesCount() const
{
    return nextCellId;
}

void SparseReconstruction::AddSinkLink(int cellId)
{
    std::map<int, std::pair<double,double> >::iterator it;
    it = sourceSinkCost.find(cellId);
    if(it != sourceSinkCost.end())
    {
        it->second.second += SINK_COST;
    }
    else
    {
        std::pair<double,double> cost;
        cost.second = SINK_COST;
        cost.first = 0;
        sourceSinkCost[cellId] = cost;
    }
}

void SparseReconstruction::GenerateGraph()
{
    maxCost = 0;
    edgeCost.clear();
    //cell2id.clear();
    nextCellId = 0;

    sourceSinkCost.clear();
    collidedCells.clear();

    //for(unsigned int i=0; i< map->vpKeyFrames.size(); i++)
    // for each point in the map
    Delaunay::Finite_vertices_iterator vhIT;
    for(vhIT = delaunayMesh.finite_vertices_begin(); vhIT != delaunayMesh.finite_vertices_end(); ++vhIT)
    {
        // this point has enough keyframes already
        Delaunay::Vertex_handle featureHandle = vhIT;
        Vector<3> featurePos;
        // is there a better way to do this?
        featurePos[0] = featureHandle->point()[0];
        featurePos[1] = featureHandle->point()[1];
        featurePos[2] = featureHandle->point()[2];

        // for each keyframe of the point
        std::vector<Vector<3> >::iterator camIT;
        for(camIT = featureHandle->info().camera_locations.begin(); camIT != featureHandle->info().camera_locations.end(); ++camIT)
        {

            // get the position of the camera
            Vector<3> CameraPos  = *camIT;

            // get the cell the camera is in
            Delaunay::Cell_handle cameraCell;
            Delaunay::Point dummy1(CameraPos[0], CameraPos[1], CameraPos[2]);
            cameraCell = delaunayMesh.locate(dummy1);


            // The cell behind the feature
            Vector<3> c2fNorm = featurePos - CameraPos;
            normalize(c2fNorm);
            Vector<3> near = featurePos + SMALLVALUE * c2fNorm;

/*            // get the cell behind the feature
            Vector<3> near = featurePos - CameraPos;
            normalize(near);
            near = featurePos + SMALLVALUE * near;
*/
            Delaunay::Point dummy2(near[0], near[1],near[2]);
            Delaunay::Cell_handle backCell = delaunayMesh.locate(dummy2);
            if(delaunayMesh.is_infinite(backCell))
            {

//                std::vector<Delaunay::Cell_handle> cells;
//                delaunayMesh.incident_cells(featureHandle,std::back_inserter(cells));

//                std::vector<Delaunay::Cell_handle> infneib;

//                for(int i=0; i<cells.size(); i++)
//                {
//                    if(delaunayMesh.is_infinite(cells[i]))
//                        infneib.push_back(cells[i]);
//                }

//                for(int i=0; i<cells.size(); i++)
//                {
//                    if(!delaunayMesh.is_infinite(cells[i]))
//                    {
//                        for(int j=0; j<infneib.size();j++)
//                        {
//                            if(cells[i]->has_neighbor(infneib[j]))
//                            {
//                                AddSinkLink(AssignIdToCell(infneib[j]));
//                                IncreaseCost(cells[i],infneib[j]);
//                            }
//                        }
//                    }
//                }

            }
            else
            {
                //if(ParamsAccess::varParams->DrawRRTTree)
                    AddSinkLink(AssignIdToCell(backCell));
            }

            Delaunay::Cell_handle curCell;
            int curId;

            // Find the cell in front of the feature
            Vector<3> f2c = -c2fNorm;

            near = featurePos + SMALLVALUE * f2c;
            Delaunay::Point dummy3(near[0], near[1],near[2]);
            curCell = delaunayMesh.locate(dummy3);

            curId = AssignIdToCell(curCell);
            Vector<3> cpoint = near;

            if(delaunayMesh.is_infinite(curCell))
            {
                AddSourceLink(curId);
                continue;
            }

            Delaunay::Cell_handle lastCell = backCell;
            int loopn = 0;
            while(true)
            {
                if(loopn++ > 1000)
                {
                    return;
                }
//                if(std::find(collidedCells.begin(),collidedCells.end(),curCell) == collidedCells.end())
//                    collidedCells.push_back(curCell);
//                else
//                {
//                    rstart = featurePos;
//                    rend = CameraPos;
//                    return;
//                }

                Delaunay::Cell_handle nextCell;
                int nextId;

                near += SMALLVALUE*f2c;

                nextCell = FindCollisionFromVertex( featurePos + SMALLVALUE * f2c, CameraPos, curCell, cpoint, lastCell);

                if(nextCell == NULL)
                    break;

                nextId = AssignIdToCell(nextCell);
                //ROS_INFO("Graph generation.. Next Id = %d", nextId);

                near = cpoint;

                IncreaseCost(nextId, curId);

                if(nextCell == cameraCell || delaunayMesh.is_infinite(nextCell))
                {
                    //ROS_INFO("-------********");
                    AddSourceLink(nextId);
                    break;
                }
                else
                {
                    //ROS_INFO("------->>>>>");
                    lastCell = curCell;
                    curCell = nextCell;
                    curId   = nextId;
                }

                //ROS_INFO("3");
            }

//            Vector<3> cameraPos = (*it)->se3CfromW.inverse().get_translation();

//            Vector<3> cp = cameraPos;
//            sf = fp;

//            Delaunay::Vertex_handle vh;
//            Delaunay::Point dummy(fp[0], fp[1], fp[2]);
//            if(!t.is_vertex(dummy,vh))
//            {
//                //ROS_INFO("??");
//                continue;
//            }

//            //start from a vertice (1) camera cell (2) a vertex on the ray to another vertex
//            Vector<3> cpoint;
//            Delaunay::Cell_handle curCell = GetCameraCell(cp,fp, cpoint);
//            if(t.is_infinite(curCell))
//                continue;

//            if(isVertex(curCell,fp))
//            {
//                //int curId;
//                //int nextId;
//                //curId = AssignIdToCell(curCell);
//                //AddSourceLink(curId);
//                continue;
//            }

//            int curId;
//            int nextId;
//            curId = AssignIdToCell(curCell);
//            AddSourceLink(curId);
//           // collidedCells.push_back(curCell);

//            Vector<3> collisionPoint = cp;

//            while(true)
//            {
//                Vector<3> near = fp-collisionPoint;
//                normalize(near);
//                near = SMALLVALUE*near + collisionPoint;

//                Delaunay::Cell_handle nextCell = FindCollisionFromVertex(near,fp,curCell,collisionPoint);
//                if(nextCell == NULL)
//                {
//                    //ROS_INFO("No Next Cell Found.");
//                    break;
//                }

//                nextId = AssignIdToCell(nextCell);
//                IncreaseCost(curId, nextId);

//                if(isVertex(nextCell,fp))
//                {
//                    Vector<3> near = fp-cameraPos;
//                    normalize(near);
//                    near = fp + SMALLVALUE * near;
//                    cpoint = near;
//                    Delaunay::Point p(near[0],near[1],near[2]);
//                    Delaunay::Cell_handle cell = t.locate(p);
//                    if(!t.is_infinite(cell))
//                    {
//                       int lastcellid = AssignIdToCell(cell);
//                       AddSinkLink(lastcellid);
//                    }
//                    else
//                    {
//                        AddSinkLink(nextId);
//                    }

//                    //ROS_INFO("At the end of the ray.");
//                    break;
//                }

//                curCell = nextCell;
//                curId = nextId;

//            }


        }

    }

//    // add all the other const weighted edges
    Delaunay::Finite_cells_iterator fci;
    for(fci = delaunayMesh.finite_cells_begin(); fci != delaunayMesh.finite_cells_end(); ++fci)
    {
        Delaunay::Cell_handle ch = fci;
        int cid = AssignIdToCell(ch);

        Delaunay::Cell_handle cn;
        for(int i=0; i<4; i++)
        {
            cn = ch->neighbor(i);

            if(delaunayMesh.is_infinite(cn))
                continue;

            int nid = AssignIdToCell(cn);

            std::pair<int,int> e;
            e.first = cid;
            e.second = nid;
            //Matrix<3,3> t;
            //Triangle3D tmp;
            //GetCommonTriangle(ch, cn, t, false, tmp.vertices, tmp.outside_vh);
            //double area = triangle_area_squared(t);
            double regFac = 0.2;//PtamParameters::varparams().MeshRegularizer /*ParamsAccess::varParams->MeshRegularizationFac*/;
            if(edgeCost.find(e) == edgeCost.end())
            {               
                edgeCost[e] = regFac;
            }
            else
            {
                edgeCost[e] += regFac;
            }

        }
    }


}

double SparseReconstruction::MinDistToVertex(Delaunay & t,  Point p)
{
    double minDist = 99999;

    if(t.number_of_vertices() < 5)
        return minDist;

    Delaunay::Vertex_handle v = t.nearest_vertex(p);
    Vector<3> pt = cgalPoint2ToonVect(p);
    if(!t.is_infinite(v))
    {
        Vector<3> vt = cgalPoint2ToonVect(v->point());
        minDist = (pt-vt)*(pt-vt);
    }

    return minDist;

}

bool SparseReconstruction::GetCommonTriangle(const Delaunay::Cell_handle &cell1, const Delaunay::Cell_handle &cell2,
                                        Matrix<3,3> &triangle, bool noCameraNeighbor,
                                        std::vector<Delaunay::Vertex_handle> &vertices, Delaunay::Vertex_handle &outside_vh)
{
    int n=0;

    for(int i=0; i<4; i++)
    {
        Delaunay::Vertex_handle v = cell1->vertex(i);
        if(cell2->has_vertex(v) && (!noCameraNeighbor || std::find(cameraVertices.begin(), cameraVertices.end(), v) == cameraVertices.end()))
        {
            vertices.push_back( v );
            Delaunay::Point p = v->point();
            triangle[n][0] = p.x();
            triangle[n][1] = p.y();
            triangle[n][2] = p.z();
            n++;
        }
        else if ( !cell2->has_vertex(v) )
        {
            outside_vh = v; //the vertex that belongs ONLY to C1 is towards the outside
        }

    }

    if(n == 3)
        return true;
    else
        return false;
}

void SparseReconstruction::DrawCell(Delaunay::Cell_handle ch)
{
    Vector<3> v[4];
    int n=0;

    if(delaunayMesh.is_infinite(ch))
    {
        return;
        for(int i=0; i<4;i++)
        {
            if(!delaunayMesh.is_infinite(ch->vertex(i)))
            {
                n++;
                copy(v[i], ch->vertex(i)->point());
            }
        }
    }
    else
    {
        for(int i=0; i<4;i++)
        {
            copy(v[i], ch->vertex(i)->point());
        }
        n = 4;
    }


    //CVD::glVertex(v[0]);
    //CVD::glVertex(v[1]);
    //CVD::glVertex(v[2]);

    if(n>3)
    {
//        CVD::glVertex(v[3]);
//        CVD::glVertex(v[1]);
//        CVD::glVertex(v[2]);

//        CVD::glVertex(v[0]);
//        CVD::glVertex(v[3]);
//        CVD::glVertex(v[2]);

//        CVD::glVertex(v[0]);
//        CVD::glVertex(v[1]);
//        CVD::glVertex(v[3]);

    }
}


void SparseReconstruction::glDraw()
{
    //ROS_INFO("check 1");

    //load the textures
//    for(int kf = 0; kf < map->vpKeyFrames.size(); kf++)
//    {
//        if(!map->vpKeyFrames[kf]->TextureLoaded)
//        {
//            //ROS_INFO("Loading Texture for KeyFrame#: %d", kf);
//            map->vpKeyFrames[kf]->TextureLoaded = true;
//            glEnable(GL_TEXTURE_2D);
//            glGenTextures(1, &map->vpKeyFrames[kf]->mnFrameTex);
//            glBindTexture(GL_TEXTURE_2D, map->vpKeyFrames[kf]->mnFrameTex);
//            glTexImage2D(GL_TEXTURE_2D,
//            0, GL_RGB,
//                            map->vpKeyFrames[kf]->imColor.size().x, map->vpKeyFrames[kf]->imColor.size().y,
//            0,GL_RGB,
//            GL_UNSIGNED_BYTE,
//            map->vpKeyFrames[kf]->imColor.data());
//            glTexParameteri(GL_TEXTURE_2D , GL_TEXTURE_MIN_FILTER , GL_NEAREST);
//            glTexParameteri(GL_TEXTURE_2D , GL_TEXTURE_MAG_FILTER , GL_NEAREST);
//        }
//    }


    if(false)
    {
        Delaunay::Finite_cells_iterator cit = delaunayMesh.finite_cells_begin();

        glLineWidth(1);

        glBegin(GL_LINES);

        for(;cit != delaunayMesh.finite_cells_end(); ++cit)
        {
            Delaunay::Vertex_handle v[4];

            v[0] = cit->vertex(0);
            v[1] = cit->vertex(1);
            v[2] = cit->vertex(2);
            v[3] = cit->vertex(3);

            Vector<3> p[4];

            copy(p[0],v[0]->point());
            copy(p[1],v[1]->point());
            copy(p[2],v[2]->point());
            copy(p[3],v[3]->point());


//            for(int i=0; i<4; i++)
//                for(int j=0; j<4; j++)
//                {
//                    CVD::glVertex(p[i]);
//                    CVD::glVertex(p[j]);
//                }
        }


        glEnd();
    }



    if(/*ParamsAccess::varParams->DrawInsideCells*/false)
    {
        for(unsigned int ii=0; ii<surfaceTriangles.size(); ii++)
        {

            double c = ((double)(ii%255))/255.0;

            glColor4f(1-c,c,1,1);
            glLineWidth(2);

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glBegin(GL_TRIANGLES);

            if(surfaceTriangles[ii]->c2 != NULL)
                DrawCell(surfaceTriangles[ii]->c2);

            glEnd();
        }

    }

    if(/*ParamsAccess::varParams->DrawSurfaceMesh*/true)
    {

        for(unsigned int ii=0; ii<surfaceTriangles.size(); ii++)
        {
            Matrix<3,3> tri = surfaceTriangles[ii]->tri;


            glColor3f(1,1,1);
            glLineWidth(3);

            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glBegin(GL_TRIANGLES);
            for(int i=0; i<3; i++)
            {
                glVertex3d(tri[i][0], tri[i][1], tri[i][2]);
            }
            glEnd();
        }

    }



//    if(/*ParamsAccess::varParams->DrawTexturedMesh*/true)
//    {
//        glColor3f(0,1,0);
//        glLineWidth(2);


//        //glPolygonMode(GL_FRONT, GL_FILL);
//        //glBegin(GL_LINES);
//        for(unsigned int ii=0; ii<surfaceTriangles.size(); ii++)
//        {
//            Matrix<3,3> tri = surfaceTriangles[ii]->tri;

//            int bestKeyframe = surfaceTriangles[ii]->KeyframeTex; //VisibilityServer::GetInstance()->GetBestKeyFrame(tri[0], tri[1], tri[2]);

//            if(bestKeyframe>=0)
//            {
//                glEnable(GL_TEXTURE_2D);
//                glBindTexture(GL_TEXTURE_2D, map->vpKeyFrames[bestKeyframe]->mnFrameTex);
//            }

//            glPolygonMode(GL_FRONT_AND_BACK , GL_FILL);

//            glBegin(GL_TRIANGLES);

//           // bool tx = TriangleNormalInCamera(map->vpKeyFrames[f[ff]]->se3CfromW, tri);
//            for(int i=0; i<3; i++)
//            {
//                glColor3f(1, 1, 1);
//                Vector<2> crd = surfaceTriangles[ii]->textureCord[i];//VisibilityServer::GetInstance()->GetImageCoord(tri[i], bestKeyframe);
//              //  if(tx)
//                {
//                    if (/*ParamsAccess::varParams->DrawTriangleValidityTrace*/false){
//                        continue;
//                        /* Here is the colour key:
//                         *      yellow triangles: too_far
//                         *      red triangles: too_big && !too_far
//                         *      green triangles: too_big && validated_by_angle && !too_far
//                         */
//                        if (surfaceTriangles[ii]->too_far)
//                        {
//                            glColor3f(0.8, 0.8, 0); //draw yellow triangle if it's too far
//                        }
//                        else
//                        {
//                            if ( surfaceTriangles[ii]->occluded )
//                            {
//                                //back facing => blue
//                                glColor3f(0, 0, 1);
//                            }
//                            else if ( surfaceTriangles[ii]->too_big )
//                            {
//                                //big & not validated => red
//                                glColor3f(1, 0, 0);
//                            }
//                            else
//                            {
//                                //not too big, not occluded and not too far => just draw the texture
//                                glTexCoord2f(crd[0], crd[1]);
//                            }
//                        }

//                    }
//                    else
//                    {
//                        if(surfaceTriangles[ii]->unexplored || bestKeyframe<0)
//                            glColor3f(0.5,0.5,0.5);
//                        else
//                            glTexCoord2f(crd[0], crd[1]);
//                    }
//                    glVertex3d(tri[i][0], tri[i][1], tri[i][2]);
//                }
//            }
//            glEnd();
//            glDisable(GL_TEXTURE_2D);


//        }

//    }

}


void SparseReconstruction::get_links_in_min_cut(int num_nodes, std::map< std::pair<int,int>, double > linksToCosts,
                                           std::map<int, std::pair<double,double> > terminal_costs,
                                           std::vector< std::pair<int,int> > &neighbor_links,
                                           std::vector< int > &source_links_in_min_cut,
                                           std::vector< int > &sink_links_in_min_cut)
{

    std::map<int, GraphType::node_id> node_id_map; //maps ids supplied as ints to "node_id"s used by the flow lib

    printf("Connecting %d nodes with %lu edges between non-terminal nodes.\n", num_nodes, linksToCosts.size());

    //populate graph
    GraphType *g = new GraphType( num_nodes, /*estimated # of edges*/ linksToCosts.size());

#ifdef _VERIFY_FLOW
    //flow verification can't be done if some edges were deleted from the supplied map
    //so we'll store them in this second structure
    std::map< std::pair<int,int>, double > revLinksToCosts;
#endif

    //add links between non-terminal nodes.
    for(std::map< std::pair<int,int>,double >::iterator mit = linksToCosts.begin(); mit != linksToCosts.end(); mit++) {

        int n1 = mit->first.first;
        int n2 = mit->first.second;

        double forward_cost = mit->second;

        //add both nodes to the graph if necessary
        GraphType::node_id n1_id = add_node_if_necessary(n1, node_id_map, g);
        GraphType::node_id n2_id = add_node_if_necessary(n2, node_id_map, g);

        //see if the link in the opposite direction exists
        std::pair<int,int> rev_pair(n2,n1);
        double reverse_cost=0.0;
        if ( linksToCosts.find(rev_pair) != linksToCosts.end() ){
            reverse_cost = linksToCosts[ rev_pair ];

            linksToCosts.erase( rev_pair );//we don't want to consider this link again, so we delete it
#ifdef _VERIFY_FLOW
            assert( revLinksToCosts.find(rev_pair)==revLinksToCosts.end() );
            revLinksToCosts[rev_pair] = reverse_cost; //only used to verify
#endif
        }

        g->add_edge(n1_id, n2_id, forward_cost, reverse_cost);
    }

    //specify costs between terminal and non-terminal nodes

    for(std::map<int, std::pair<double,double> >::iterator terms_it = terminal_costs.begin(); terms_it != terminal_costs.end(); terms_it++){

        //nodes that have no non-terminal neighbours will have to be added
        GraphType::node_id id = add_node_if_necessary(terms_it->first, node_id_map, g);

        double costFromSource = terms_it->second.first;
        double costToSink = terms_it->second.second;

        g -> add_tweights( id, costFromSource, costToSink );
    }

    //calc flow
    double flow = g -> maxflow();

 #ifdef _VERIFY_FLOW
    double min_cut_capacity=0.0;
 #endif


    //we now find which edges have one end in the sink set and another in the source set
    for(std::map< std::pair<int,int>,double >::iterator mit = linksToCosts.begin(); mit != linksToCosts.end(); mit++)
    {
        int n1 = mit->first.first;
        int n2 = mit->first.second;

        //get IDs
        GraphType::node_id n1_id = node_id_map[n1];
        GraphType::node_id n2_id = node_id_map[n2];

#ifdef _VERIFY_FLOW
        std::pair<int,int> rev_pair(n2, n1);
        if (g->what_segment(n1_id)==GraphType::SINK   && g->what_segment(n2_id)==GraphType::SOURCE){
            //neighbor_links.push_back( mit->first ); //this edge IS NOT in the min cut!!!! (but the reverse edge is, if it exists)

            //the min cut considers only edges from the source set to the sink set
            if ( revLinksToCosts.find( rev_pair ) != revLinksToCosts.end() ){
                min_cut_capacity += revLinksToCosts[rev_pair];
                neighbor_links.push_back( rev_pair );
            }

        } else if (g->what_segment(n1_id)==GraphType::SOURCE && g->what_segment(n2_id)==GraphType::SINK) {
            neighbor_links.push_back( mit->first );

            min_cut_capacity += mit->second;
        }
#else
        if (   (g->what_segment(n1_id)==GraphType::SINK   && g->what_segment(n2_id)==GraphType::SOURCE)
            || (g->what_segment(n1_id)==GraphType::SOURCE && g->what_segment(n2_id)==GraphType::SINK) )
        {
            neighbor_links.push_back( mit->first );
        }
 #endif
    }

    //put nodes that have links to the terminal nodes
    for(std::map<int, std::pair<double,double> >::iterator t_it = terminal_costs.begin(); t_it != terminal_costs.end(); t_it++)
    {
        int n = t_it->first;
        GraphType::node_id id = node_id_map[n];

        double costFromSource = t_it->second.first;
        double costToSink = t_it->second.second;

        /*          2
         * (source)---->(id)
         * If the (id) node is in the sink set, then this link is in the cut
         */

        if (g->what_segment(id)==GraphType::SINK && costFromSource>0)
        {
            source_links_in_min_cut.push_back(n);
#ifdef _VERIFY_FLOW
            min_cut_capacity += costFromSource;
#endif
        }
        else if ( g->what_segment(id)==GraphType::SOURCE && costToSink>0)
        {
            sink_links_in_min_cut.push_back(n);
#ifdef _VERIFY_FLOW
            min_cut_capacity += costToSink;
#endif
        }
    }

#ifdef _VERIFY_FLOW

    printf("There are %lu links with one end in the sink set and the other in the source set, "
           "\n\t %lu links connected to the source in the min cut, "
           "\n\t %lu links connected to the sink in the min cut\n", neighbor_links.size(), source_links_in_min_cut.size(), sink_links_in_min_cut.size());
    printf("Calculated flow to be %f and summed edges in min cut to be %f.\n", flow, min_cut_capacity);
    assert(fabs(flow-min_cut_capacity)<0.0001);
#endif

    delete g;
}

SparseReconstruction::GraphType::node_id SparseReconstruction::add_node_if_necessary(int n, std::map<int, GraphType::node_id> &node_id_map, GraphType *g)
{
    GraphType::node_id  id;
    if ( node_id_map.find(n) != node_id_map.end() )
    {
        id = node_id_map[ n ];
    }
    else
    {
        id = g->add_node();
        node_id_map[ n ] = id;
    }
    return id;
}

template<typename VertexListGraph, typename Vertex>
std::map<Vertex, SparseReconstruction::Triangle3D> SparseReconstruction::map_vertices_to_triangles(std::vector<Triangle3D> tris, VertexListGraph &g)
{
    using namespace boost;
    using namespace std;

    typedef typename graph_traits<VertexListGraph>::vertex_iterator VItr;

    std::map<Vertex, Triangle3D> v_pmap;

    VItr vi, ve;
    int idx(-1);
    for (boost::tie(vi, ve) = vertices(g); vi != ve; ++vi)
    {
        Vertex v(*vi);
        if (idx != -1) //we leave the first vertex unmapped, since it's a dummy vertex
        {
            v_pmap[v] = tris[idx];
        }

        idx++;
    }

    return v_pmap;
}

template<typename WeightMap, typename Vertex, typename VertexListGraph>
void SparseReconstruction::create_connected_graph(VertexListGraph &g, WeightMap wmap, std::map<Vertex,Triangle3D> vmap)
{
    using namespace std;
    using namespace boost;
    typedef typename graph_traits<VertexListGraph>::edge_descriptor Edge;
    typedef typename graph_traits<VertexListGraph>::vertex_iterator VItr;

    Edge e;
    bool inserted;

    pair<VItr, VItr> verts(vertices(g));
    for (VItr src(verts.first); src != verts.second; src++)
    {
        for (VItr dest(src); dest != verts.second; dest++)
        {
            if (dest != src)
            {
                double weight = 0.0;

                //the weights to and from the dummy vertex (indexes NULL) will be left as zero
                if (vmap.find(*src) != vmap.end() && vmap.find(*dest) != vmap.end()){
                    Vector<3> src_center = get_triangle_centroid( vmap[*src].tri );
                    Vector<3> dest_center = get_triangle_centroid( vmap[*dest].tri );

                    weight = euclidean_distance_squared(src_center, dest_center);
                }

                boost::tie(e, inserted) = add_edge(*src, *dest, g);

                wmap[e] = weight;
            }


        }
    }
    ROS_INFO("Found %lu vertices and %lu edges\n", num_vertices(g), num_edges(g));
}


SparseReconstruction::CellInfo::CellInfo() : id(CellInfo::NO_ID) {}

SparseReconstruction::Triangle3D* SparseReconstruction::CellInfo::getCommonTriangle3D( int cell_id )
{
    for (unsigned int i=0; i < surfaceTris.size(); i++)
    {
        if (surfaceTris[i]->c1->info().id == cell_id
                || surfaceTris[i]->c2->info().id == cell_id )
            return surfaceTris[i];
    }
    return NULL;
}

unsigned long SparseReconstruction::VertexInfo::nextVertexId(0);

SparseReconstruction::VertexInfo::VertexInfo()
{
    id = nextVertexId++;
}



//bool SparseReconstruction::saveMesh(const std::string &dir)
//{
//    std::string dirName = dir;
//    if(*dirName.rbegin() != '/')
//        dirName += '/';
//    std::set<int> framesUsed;
//    //char curDIR[256];
//    //getcwd(&curDIR, 256);
//    std::string xmlFile = dirName + "mapdata.xml";

//    ROS_INFO("Saving map to %s", dirName.c_str());

//    struct stat dirStatus;
//    if(!(stat(dirName.c_str(),&dirStatus) == 0 && (dirStatus.st_mode & S_IFDIR) != 0))
//    {
//        ROS_INFO("Directory '%s' doesn't exist.  Creating...", dirName.c_str());
//        // doesn't exist, make it
//        if(mkdir(dirName.c_str(), S_IRWXU | S_IRWXG | S_IRWXO ) != 0)
//        {
//            ROS_WARN("Failed to create %s! Aborting!", dirName.c_str());
//            return false;
//        }
//    }else
//    {
//        ROS_WARN("Directory already '%s' exists!  Aborting!", dirName.c_str());
//        return false;
//        struct dirent *fileHandle;
//        DIR *dirHandle;

//        dirHandle = opendir(dirName.c_str());

//        while ( (fileHandle = readdir(dirHandle)) )
//        {
//            if(fileHandle->d_type != DT_DIR)
//            {
//                // build the full path for each file in the folder
//                std::string fName(fileHandle->d_name);
//                fName = dirName + fName;
//                //ROS_INFO("pretending to remove %s", fName.c_str());
//                remove(fName.c_str());
//            }
//        }
//    }

//    std::ofstream xmlf(xmlFile.c_str(), std::ofstream::binary | std::ofstream::trunc);
//    if(!xmlf || !xmlf.is_open())
//        return false;

//    ROS_INFO("Writing XML...");
//    // NOTE: faking XML as there is little reason to use an xml library if this is only going to be replaced eventually

//    /* format:
//     *      <map>
//     *          <version='0.0.0'/>
//     *          <triangles>
//     *              ...
//     *          </triangles>
//     *          <frames>
//     *              ...
//     *          </frames>
//     *      </map>*/

//    // writing header

//    xmlf << "<map>\n";
//    xmlf << "\t<version value='0.0.0'/>\n";
//    xmlf << "\t<triangles count='" << surfaceTriangles.size() <<"'>\n";

//    // then for each triangle in the mesh
//    for(unsigned int ii=0; ii<surfaceTriangles.size(); ii++)
//    {
//        // get the vertexes of the triangle
//        Matrix<3,3> tri = surfaceTriangles[ii]->tri;

//        // and the texture id for the texture to use
//        int bestKeyFrame = surfaceTriangles[ii]->KeyframeTex;

//        // add to used frames
//        framesUsed.insert(bestKeyFrame);

//        // write triangle out to file
//        /* format:
//         *          <tri frame = bestKeyFrame>
//         *              <vertex>
//         *                  <pos x='x' y='y' z='z'/>
//         *                  <texture x='x' y='y'/>
//         *              </vertex>
//         *              <vertex> ... x2
//         *              </vertex>
//         *          </tri>*/

//        xmlf << "\t\t<tri frame='" << ((bestKeyFrame<0)?0:bestKeyFrame) << "'>\n";

//        // then for each vertex in triangle
//        for(unsigned int i=0; i<3; i++)
//        {
//            xmlf << "\t\t\t<vertex>\n";
//            // get the texture coordinate
//            Vector<2> crd = surfaceTriangles[ii]->textureCord[i];
//            // and the vertex pos
//            Vector<3> vPos = tri[i];

//            xmlf << "\t\t\t\t<pos x='" << vPos[0] << "' y='" << vPos[1] << "' z='" << vPos[2] << "'/>\n";
//            xmlf << "\t\t\t\t<texture x='" << crd[0] << "' y='" << crd[1] << "'/>\n";

//            xmlf << "\t\t\t</vertex>\n";
//        }

//        xmlf << "\t\t</tri>\n";
//    }

//    // closing triangles
//    xmlf << "\t</triangles>\n";
//    // openning frames
//    xmlf << "\t<frames count='"<< framesUsed.size() <<"'>\n";

//    ROS_INFO("Writing keyframes...");

//    // Serialize textures in file
//    for(std::set<int>::iterator I = framesUsed.begin(); I != framesUsed.end(); ++I)
//    {
//        int kf = *I;
//        if(kf < 0) continue;

//        std::stringstream ss;
//        ss << kf;
//        std::string frameName = "frame" + ss.str() + ".jpg";
//        std::string frameFile = dirName + frameName;
//        // going to save images outside of the xml with ( http://www.edwardrosten.com/cvd/cvd/html/group__gImageIO.html#ga7b74b54ba166c9f3901e5eec9fa07bfe )
//        // then link the xml to the image name for loading purposes

//        std::ofstream frameF(frameFile.c_str(), std::ofstream::binary | std::ofstream::trunc);
//        if(!frameF || !frameF.is_open())
//        {
//            // close out and fail
//            xmlf.close();
//            return false;
//        }
//        // otherwise write the frame to disk
//        // TODO, write to disk

//        CVD::img_save(map->vpKeyFrames[kf]->imColor, frameF, CVD::ImageType::JPEG );


//        // done writing image
//        frameF.close();
//        // write frame info to xml.
//        /* format:
//         *      <frame id ='kf' file='file'/>
//         */

//        xmlf << "\t\t<frame id='" << kf << "' file='" << frameName << "'/>\n";
//    }

//    // closing textures and opening vertexs
//    xmlf << "\t</frames>\n";
//    // and closing root node
//    xmlf << "</map>\n";

//    xmlf.close();
//    return true;
//}

//bool SparseReconstruction::savePoints(const std::string &file)
//{
//    std::string xmlFile = file + ".xml";

//    ROS_INFO("Saving points to %s", xmlFile.c_str());

//    std::ofstream xmlf(xmlFile.c_str(), std::ofstream::binary | std::ofstream::trunc);
//    if(!xmlf || !xmlf.is_open())
//        return false;

//    long unsigned int keyFramePointSize = 0;

//    ROS_INFO("Writing XML...");
//    // NOTE: faking XML as there is little reason to use an xml library if this is only going to be replaced eventually

//    /* format:
//     *      <pointmap>
//     *          <version='0.0.0'/>
//     *          <points count=..>
//     *              <point x=.. y=.. z=..>
//     *                  <kfs count=..>
//     *                      <kf id=../>
//     *                      ...
//     *                  </kfs>
//     *              </point>
//     *              ...
//     *          </points>
//     *          <keyframes count=..>
//     *              <kf id=.. x=.. y=.. z=.. />
//     *          </keyframes>
//     *      </map>*/

//    // writing header

//    xmlf << std::fixed << std::setprecision(10);

//    xmlf << "<pointmap>\n";
//    xmlf << "\t<version value='0.0.0'/>\n";
//    xmlf << "\t<points count='" << map->vpPoints.size() <<"'>\n";

//    std::map<KeyFrame::Ptr, unsigned int> kf2id;
//    unsigned int kfidcount = 0;

//    // then for each point in the map
//    for(unsigned int ii=0; ii < map->vpPoints.size(); ii++)
//    {
//        boost::shared_ptr<MapPoint> cur = map->vpPoints[ii];
//        // get the position of the point
//        Vector<3> pos = cur->v3WorldPos;

//        xmlf << "\t\t<point x='" << pos[0] << "' y='" << pos[1] << "' z='" << pos[2] << "'>\n";

//        // get the number of keyframes that saw it
//        unsigned int kfc = cur->pMMData->sMeasurementKFs.size();

//        xmlf << "\t\t\t<kfs count='" << kfc << "'>\n";

//        // for each keyframe
//        for(std::set<KeyFrame::Ptr>::iterator iitt = cur->pMMData->sMeasurementKFs.begin(); iitt != cur->pMMData->sMeasurementKFs.end(); ++iitt)
//        {
//            boost::shared_ptr<KeyFrame> curkf = *iitt;
//            if(kf2id.find(curkf) == kf2id.end())
//            {
//                // assign a new id
//                kf2id[curkf] = ++kfidcount;
//            }

//            xmlf << "\t\t\t\t<kf id='" << kf2id[curkf] << "'/>\n";
//            keyFramePointSize++;
//        }

//        // closing keyframes
//        xmlf << "\t\t\t</kfs>\n";

//        // closing point
//        xmlf << "\t\t</point>\n";
//    }

//    ROS_INFO("Saved %lu points", map->vpPoints.size());

//    ROS_INFO("Saved %lu unique keyframe-point pairs", keyFramePointSize);

//    // closing points
//    xmlf << "\t</points>\n";

//    xmlf << "\t<keyframes count='" << kfidcount << "'>\n";

//    // for each keyframe
//    for(std::map<KeyFrame::Ptr, unsigned int>::iterator iitt = kf2id.begin(); iitt != kf2id.end(); ++iitt)
//    {
//        KeyFrame::Ptr curkf = (*iitt).first;
//        unsigned int kfid = (*iitt).second;
//        Vector<3> kfpos  = curkf->se3CfromW.inverse().get_translation();

//        xmlf << "\t\t<kf id='" << kfid << "' x='" << kfpos[0] << "' y='" <<kfpos[1] << "' z='" << kfpos[2] << "'/>\n";
//    }

//    ROS_INFO("Saved %lu keyframes", kf2id.size());

//    // closing keyframes
//    xmlf << "\t</keyframes>\n";

//    // closing map
//    xmlf << "</pointmap>\n";

//    long unsigned int fileSize = static_cast<unsigned int>(xmlf.tellp());
//    ROS_INFO("Saved %lu bytes", fileSize);

//    xmlf.close();
//    return true;
//}
