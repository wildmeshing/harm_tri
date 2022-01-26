#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Real_timer.h>

#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_ply_points.h>

#include "CLI11.hpp"

#include <iostream>
#include <fstream>
#include <vector>
//#include <set>
//#include <unordered_set>
//#include <bitset>
#include <queue>
//#include <map>
//#include <chrono>
//#include <random>
#include <algorithm>

//#include <boost/math/tools/minima.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;

typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb;
typedef CGAL::Triangulation_cell_base_with_info_3<int, K>       Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;

typedef CGAL::Delaunay_triangulation_3<K,Tds>      DT;
typedef CGAL::Triangulation_3<K,Tds>      Triangulation;

typedef Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
typedef Triangulation::Finite_edges_iterator Finite_edges_iterator;
typedef Triangulation::Finite_facets_iterator Finite_facets_iterator;
typedef Triangulation::Finite_cells_iterator Finite_cells_iterator;

typedef Triangulation::Point          Point;
typedef Triangulation::Facet          Facet;
typedef Triangulation::Edge          Edge;
typedef Triangulation::Vertex_handle          VH;
typedef Triangulation::Cell_handle          CH;
typedef Triangulation::Cell_circulator      CC;
typedef Triangulation::Facet_circulator      FC;


//typedef K::Vector_3 Vector;

struct flip_info {
    CH c; // the cell defining the edge or face to be flipped
    int i,j; // cell indices defining edge or face (j=-1)
    int ci1, ci2, ci3; // cell indices (info) at time of creation
    FT hi1,hi2,hi3; // the 2 or 3 harmonic indices for the cells after flipping
                    // so we don't have to recompute them
    FT val; // value for priority queue
};

struct flip_comp {
    bool operator()(flip_info a,flip_info b) const {
        return a.val < b.val;
    }
};

std::vector<FT> hi;


int n, cci = 0;
Triangulation T;

//std::vector<bool> bdy;
//
//Polyhedron bdy_poly;
//Tree bdy_tree;
//bool fix_bdy_points;
//
//static int histcount = 0;


void write_cgal(std::string filename)
{
    std::cout << "Writing tets to " << filename << std::endl;
    std::ofstream of(filename);
    of << T;
    of.close();
}

void write_tetgen(std::string filename)
{
    // tetget/tetview really expects the points to be sorted by index...
    std::vector<VH> V(n);
    Finite_vertices_iterator vit = T.finite_vertices_begin();
    for(; vit != T.finite_vertices_end(); ++vit)
        V[vit->info()]=vit;

    std::ofstream node(filename+".node");
    node << n << " 3 0 0" << std::endl;
    for (int i = 0; i < n; i++)
    {
        node << (i+1);
        Point p = V[i]->point();
        for (int j = 0; j < 3; j++)
            node << " " << p[j];
        node << std::endl;
    }

    node.close();
    std::ofstream ele(filename+".ele");
    ele << T.number_of_finite_cells() << " 4 0" << std::endl;
    Finite_cells_iterator cit = T.finite_cells_begin();
    int cid = 0;
    for (; cit != T.finite_cells_end(); ++cit)
    {
        ele << (++cid);
        for (int j = 0; j < 4; j++)
            ele << " " << (cit->vertex(j)->info()+1);
        ele << std::endl;
    }
    ele.close();
}


void stats()
{
    std::vector<double> da,vol,csr,isr,ar;
    double hti = 0.0;
    
    Finite_cells_iterator cit = T.finite_cells_begin();
    int nc = 0;
    for (; cit != T.finite_cells_end(); ++cit)
    {
        std::array<CGAL::Vector_3<K>,4> normal;
        double a = 0.0, a2 = 0.0;
        
        for (int i = 0; i < 4; i++)
        {
            normal[i] = CGAL::unit_normal(T.point(cit,i),
                                          T.point(cit,(i+1)%4),
                                          T.point(cit,(i+2)%4));
            if (i&1) normal[i] *= -1;
            double fa2 = CGAL::squared_area(T.point(cit,i),
                                            T.point(cit,(i+1)%4),
                                            T.point(cit,(i+2)%4));
            a += sqrt(fa2);
            a2 += fa2;
 
        }
        double v = T.tetrahedron(cit).volume();
        if (v < 0) std::cerr << "Tet with negative volume: " << cit->vertex(0)->info() << std::endl;
        vol.push_back(v);
        hti += a2/v;
        double r = 3.0*v/a;
        isr.push_back(r);

        double R = sqrt(CGAL::squared_radius(T.point(cit,0),
                                             T.point(cit,1),
                                             T.point(cit,2),
                                             T.point(cit,3)));
        csr.push_back(R);
        
        ar.push_back(R/r);

        for (int i = 0; i < 3; i++)
            for (int j = i+1; j < 4; j++)
                da.push_back(acos(-(normal[i]*normal[j]))*180.0/M_PI);
        nc++;
    }
    std::sort(da.begin(), da.end());
    std::cout << "Number of vertices: " << T.number_of_vertices() << std::endl;
    std::cout << "Number of finite cells: " << nc << std::endl;
    std::cout << "Trace of Laplacian: " << (hti) << std::endl;
    std::cout << "Mean Harmonic index: " << (hti / double(nc)) << std::endl;
    std::cout << std::endl << "Dihedral angles" << std::endl;
    std::cout << "Min: " << da.front() << std::endl;
    std::cout << "5th percentile: " << da[da.size()/20] << std::endl;
    std::cout << "Median: " << da[da.size()/2] << std::endl;
    std::cout << "95th percentile: " << da[19*da.size()/20] << std::endl;
    std::cout << "Max: " << da.back() << std::endl;

    std::sort(ar.begin(),ar.end());
    std::cout << std::endl << "Aspect ratio (circumradius / inradius)" << std::endl;
    std::cout << "5th percentile: " << ar[ar.size()/20] << std::endl;
    std::cout << "Median: " << ar[ar.size()/2] << std::endl;
    std::cout << "95th percentile: " << ar[19*ar.size()/20] << std::endl;

    /*
    double mir = 0.0;
    for (auto ir : isr) mir += ir;
    std::cerr << std::endl << "Mean inradius: " << (mir/double(isr.size())) << std::endl;

    double mcr = 0.0;
    for (auto cr : csr) mcr += cr;
    std::cerr << std::endl << "Mean circumradius: " << (mcr/double(csr.size())) << std::endl;
    */
}


// the following two functions are taken from CGAL source because they were private
void set_adjacency(CH c0, int i0,
                   CH c1, int i1)
{
    c0->set_neighbor(i0,c1);
    c1->set_neighbor(i1,c0);
}


void change_orientation(CH c)
{
    VH tmp_v = c->vertex(0);
    c->set_vertex(0, c->vertex(1));
    c->set_vertex(1, tmp_v);
    CH tmp_c = c->neighbor(0);
    c->set_neighbor(0, c->neighbor(1));
    c->set_neighbor(1, tmp_c);
}



FT eta(CH ch)
{
    FT sas(0);
    for (int i = 0; i < 4; i++)
        sas += T.triangle(ch,i).squared_area();
    return sas / T.tetrahedron(ch).volume();
}


bool check_face(CH c, int i, flip_info &fi)
{
    CH n = c->neighbor(i);
    // get the 2 vertex handles opposite the face, stop if any is infinite
    VH vhl = c->vertex(i);
    if (T.is_infinite(vhl)) return false;
    int in = n->index(c);
    VH vhr = n->vertex(in);
    if (T.is_infinite(vhr)) return false;
    
    // get the other three vertices and their handles
    int i1 = (i+1)&3;
    int i2 = (i+2)&3;
    int i3 = (i+3)&3;
    
    VH vhf1 = c->vertex(i1);
    VH vhf2 = c->vertex(i2);
    VH vhf3 = c->vertex(i3);
    
    // the three new tets are
    // vhl, vhf1, vhf2, vhr
    // vhl, vhf2, vhf3, vhr
    // vhl, vhf3, vhf1, vhr
    
    // compute their volumes, if any is neagtive stop
    
    FT vol1 = CGAL::volume(vhl->point(),
                              vhf1->point(),
                              vhf2->point(),
                              vhr->point());
    if (vol1 <= 0.0) return false;
    
    FT vol2 = CGAL::volume(vhl->point(),
                              vhf2->point(),
                              vhf3->point(),
                              vhr->point());
    if (vol2 <= 0.0) return false;
    
    FT vol3 = CGAL::volume(vhl->point(),
                              vhf3->point(),
                              vhf1->point(),
                              vhr->point());
    if (vol3 <= 0.0) return false;
    
    // compute the squared face areas
    
    FT sasl1 = CGAL::squared_area(vhl->point(),vhf1->point(),vhf2->point());
    FT sasl2 = CGAL::squared_area(vhl->point(),vhf2->point(),vhf3->point());
    FT sasl3 = CGAL::squared_area(vhl->point(),vhf3->point(),vhf1->point());
    
    FT sasr1 = CGAL::squared_area(vhr->point(),vhf1->point(),vhf2->point());
    FT sasr2 = CGAL::squared_area(vhr->point(),vhf2->point(),vhf3->point());
    FT sasr3 = CGAL::squared_area(vhr->point(),vhf3->point(),vhf1->point());
    
    FT sasi1 = CGAL::squared_area(vhl->point(),vhf1->point(),vhr->point());
    FT sasi2 = CGAL::squared_area(vhl->point(),vhf2->point(),vhr->point());
    FT sasi3 = CGAL::squared_area(vhl->point(),vhf3->point(),vhr->point());
    
    fi.hi1 = (sasl1+sasr1+sasi1+sasi2)/vol1;
    fi.hi2 = (sasl2+sasr2+sasi2+sasi3)/vol2;
    fi.hi3 = (sasl3+sasr3+sasi3+sasi1)/vol3;
    
    fi.val = hi[c->info()]+hi[n->info()] - (fi.hi1+fi.hi2+fi.hi3);
    if (fi.val <= FT(0)) return false;

    fi.c = c;
    fi.i = i;
    fi.j = -1;
    fi.ci1 = c->info();
    fi.ci2 = n->info();
    fi.ci3 = -1;
    return true;
}

bool perform_face_flip(flip_info fi, Edge &e)
{
    CH c = fi.c;
    int i = fi.i;
    if (c==NULL || c->info() != fi.ci1) return false;
    CH n = c->neighbor(i);
    if (n==NULL || n->info() != fi.ci2) return false;

    // get the 2 vertex handles opposite the face
    VH vhl = fi.c->vertex(i);
    int in = n->index(c);
    VH vhr = n->vertex(in);
    
    // get the other three vertices
    int i1 = (i+1)&3;
    int i2 = (i+2)&3;
    int i3 = (i+3)&3;
    
    int in1 = n->index(c->vertex(i1));
    int in2 = n->index(c->vertex(i2));
    int in3 = n->index(c->vertex(i3));
    
    T.tds().set_adjacency(c, i, n->neighbor(in3), n->neighbor(in3)->index(n));
    c->set_vertex( i3, n->vertex(in) );
    
    T.tds().set_adjacency(n, in, c->neighbor(i1), c->neighbor(i1)->index(c));
    n->set_vertex( in1, c->vertex(i) );
    
    CH cnew = T.tds().create_cell(c->vertex(i), c->vertex(i1),
                            n->vertex(in), n->vertex(in3));
    
    T.tds().set_adjacency(cnew, 0, n->neighbor(in2), n->neighbor(in2)->index(n));
    T.tds().set_adjacency(cnew, 1, n, in2);
    T.tds().set_adjacency(cnew, 2, c->neighbor(i2), c->neighbor(i2)->index(c));
    T.tds().set_adjacency(cnew, 3, c, i2);
    T.tds().set_adjacency(c, i1, n, in3);
    
    if ((i&1) != 0)
        change_orientation(cnew);
    
    c->vertex(i1)->set_cell(cnew);
    c->vertex(i2)->set_cell(c);
    n->vertex(in3)->set_cell(n);

    // cleanup

    c->info() = cci++;
    hi.push_back(fi.hi1);
    
    n->info() = cci++;
    hi.push_back(fi.hi2);

    cnew->info() = cci++;
    hi.push_back(fi.hi3);

    e.first = c;
    e.second = i;
    e.third = i3;
    
    return true;
}


bool check_edge(CH c, int i, int j, flip_info &fi)
{
    // get vertex handles of vertices around edge
    int next = CGAL::Triangulation_utils_3::next_around_edge(i,j);

    CH c1 = c->neighbor( next );
    VH vhf1 = c->vertex( next ); // will become vertex of c1
    if (T.is_infinite(vhf1)) return false;
    
    int i1 = c1->index( c->vertex(i) );
    int j1 = c1->index( c->vertex(j) );
    int next1 = CGAL::Triangulation_utils_3::next_around_edge(i1,j1);
    
    CH c2  = c1->neighbor( next1 );
    VH vhf2 = c1->vertex( next1 ); // will become vertex of c3
    if (T.is_infinite(vhf2)) return false;

    int i2 = c2->index( c->vertex(i) );
    int j2 = c2->index( c->vertex(j) );
    int next2 = CGAL::Triangulation_utils_3::next_around_edge(i2,j2);
    
    if (c != c2->neighbor( next2 ) ) return false; // edge degree != 3
    VH vhf3 = c2->vertex( next2 );
    if (T.is_infinite(vhf3)) return false;

    VH vhl = c->vertex(i);
    VH vhr = c->vertex(j);
    
    // the two new tets are
    // vhl, vhf1, vhf2, vhf3
    // vhf1, vhf2, vhf3, vhr

    // compute their volumes, if any is neagtive stop
    
    FT voll = CGAL::volume(vhl->point(),
                              vhf1->point(),
                              vhf2->point(),
                              vhf3->point());
    if (voll <= 0.0) return false;
    
    FT volr = CGAL::volume(vhf1->point(),
                              vhf2->point(),
                              vhf3->point(),
                              vhr->point());
    if (volr <= 0.0) return false;
    
    // compute the squared face areas
    
    FT sasl1 = CGAL::squared_area(vhl->point(),vhf1->point(),vhf2->point());
    FT sasl2 = CGAL::squared_area(vhl->point(),vhf2->point(),vhf3->point());
    FT sasl3 = CGAL::squared_area(vhl->point(),vhf3->point(),vhf1->point());
    
    FT sasr1 = CGAL::squared_area(vhr->point(),vhf1->point(),vhf2->point());
    FT sasr2 = CGAL::squared_area(vhr->point(),vhf2->point(),vhf3->point());
    FT sasr3 = CGAL::squared_area(vhr->point(),vhf3->point(),vhf1->point());
    
    FT sasi = CGAL::squared_area(vhf1->point(),vhf2->point(),vhf3->point());
    
    fi.hi1 = (sasl1+sasl2+sasl3+sasi)/voll;
    fi.hi2 = (sasr1+sasr2+sasr3+sasi)/volr;
    fi.hi3 = FT(0);
    
    if (hi[c1->info()]+hi[c2->info()]+hi[c->info()] <= (fi.hi1+fi.hi2))
        return false;

    fi.val = hi[c1->info()]+hi[c2->info()]+hi[c->info()] - (fi.hi1+fi.hi2);
    fi.c = c;
    fi.i = i;
    fi.j = j;
    fi.ci1 = c->info();
    fi.ci2 = c1->info();
    fi.ci3 = c2->info();
    return true;
}

bool perform_edge_flip(flip_info fi, Facet &f)
{
    CH c = fi.c;
    int i = fi.i, j = fi.j;
    
    if (c==NULL || c->info() != fi.ci1) return false;
    // get vertex handles of vertices around edge
    int next = CGAL::Triangulation_utils_3::next_around_edge(i,j);
    
    CH c1 = c->neighbor( next );
    if (c1==NULL || c1->info() != fi.ci2) return false;
    VH v1 = c->vertex( next ); // will become vertex of c1
    
    int i1 = c1->index( c->vertex(i) );
    int j1 = c1->index( c->vertex(j) );
    int next1 = CGAL::Triangulation_utils_3::next_around_edge(i1,j1);
    
    CH c2  = c1->neighbor( next1 );
    if (c2==NULL || c2->info() != fi.ci3) return false;
    VH v2 = c1->vertex( next1 ); // will become vertex of c2
    
    int i2 = c2->index( c->vertex(i) );
    int j2 = c2->index( c->vertex(j) );
    int next2 = CGAL::Triangulation_utils_3::next_around_edge(i2,j2);
    
    VH v3 = c2->vertex( next2 );
    
    c->vertex(i)->set_cell(c1);
    c->vertex(j)->set_cell(c2);
    
    c1->set_vertex( j1, v1 );
    v1->set_cell(c1);
    c2->set_vertex( i2, v2 );
    v2->set_cell(c2);
    
    set_adjacency(c1, next1,c2->neighbor(j2), c2->neighbor(j2)->index(c2));
    set_adjacency(c2,c2->index(v1),c1->neighbor(i1),c1->neighbor(i1)->index(c1));
    
    set_adjacency(c1, i1, c2, j2);
    
    set_adjacency(c1, 6-i1-j1-next1, c->neighbor(j), c->neighbor(j)->index(c));
    set_adjacency(c2, next2, c->neighbor(i), c->neighbor(i)->index(c));
    
    v3->set_cell( c2 );
    
    // cleanup
    
    T.tds().delete_cell( c ); // perhaps better not to do this? or is it necessary?
    c->info() = -1;
    
    c1->info() = cci++;
    hi.push_back(fi.hi1);
    
    c2->info() = cci++;
    hi.push_back(fi.hi2);
    
    f.first = c1;
    f.second = i1;
    
    return true;
}



void flip2harmonic(bool dofaces)
{
    CGAL::Real_timer rt;
    double total_time = 0.0;
    
    hi.resize(0);
    
    // (re)compute harmonic indices, set up original cell ids
    
    rt.start();
    Finite_cells_iterator cit = T.finite_cells_begin();
    for (cci = 0; cit != T.finite_cells_end(); ++cit, ++cci)
    {
        cit->info() = cci;
        hi.push_back(eta(cit));
    }
    rt.stop();
    std::cerr << "Computed " << cci << " harmonic indices in " << rt.time() << " seconds" << std::endl;
    total_time += rt.time();
    
    std::priority_queue<flip_info, std::vector<flip_info>, flip_comp> fpq;
    int nf = 0, ne = 0;
    
    // check all finite edges, insert the ones that should be flipped into queue
    rt.reset(); rt.start();
    Finite_edges_iterator eit = T.finite_edges_begin();
    for (; eit != T.finite_edges_end(); ++eit)
    {
        flip_info fi;
        if (check_edge(eit->first,eit->second,eit->third,fi))
        {
            fpq.push(fi);
            ne++;
        }
    }
    rt.stop();
    std::cerr << "Found " << ne << " initial edge flips in " << rt.time() << " seconds" << std::endl;
    total_time += rt.time();

    // check all finite faces, insert the ones that should be flipped into queue
    if (dofaces)
    {
        rt.reset(); rt.start();
        Finite_facets_iterator fit = T.finite_facets_begin();
        for (; fit != T.finite_facets_end(); ++fit)
        {
            flip_info fi;
            if (check_face(fit->first,fit->second,fi))
            {
                fpq.push(fi);
                nf++;
            }
        }
        rt.stop();
        std::cerr << "Found " << nf << " initial face flips in " << rt.time() << " seconds" << std::endl;
        total_time += rt.time();
    }
    
    rt.reset(); rt.start();
    nf = 0; ne = 0;
    while (!fpq.empty())
    {
        flip_info fi = fpq.top();
        fpq.pop();
        
        // if j entry in face_info is -1 this is a face flip
        if (fi.j < 0)
        {
            Edge e;
            if (perform_face_flip(fi,e))
            {
                // go through three tets incident on new edge
                // for each of them check faces on convex hull
                // and the one edge opposite to the new diagonal
                
                CH c = e.first;
                int i = e.second;
                int j = e.third;
                int next = CGAL::Triangulation_utils_3::next_around_edge(i,j);
                int k = 6-(i+j+next);
 
                flip_info fifi, fifj, fie;
                if (check_face(c,i,fifi))
                    fpq.push(fifi);
                if (check_face(c,j,fifj))
                    fpq.push(fifj);
                if (check_edge(c,next,k,fie))
                    fpq.push(fie);
                
                
                CH c1 = c->neighbor( next );
                VH v1 = c->vertex( next );
                
                int i1 = c1->index( c->vertex(i) );
                int j1 = c1->index( c->vertex(j) );
                int next1 = CGAL::Triangulation_utils_3::next_around_edge(i1,j1);
                int k1 = 6-(i1+j1+next1);

                flip_info fifi1, fifj1, fie1;
                if (check_face(c,i,fifi1))
                    fpq.push(fifi1);
                if (check_face(c,j,fifj1))
                    fpq.push(fifj1);
                if (check_edge(c,next,k,fie1))
                    fpq.push(fie1);
                
                CH c2  = c1->neighbor( next1 );
                VH v2 = c1->vertex( next1 );
                
                int i2 = c2->index( c->vertex(i) );
                int j2 = c2->index( c->vertex(j) );
                int next2 = CGAL::Triangulation_utils_3::next_around_edge(i2,j2);
                int k2 = 6-(i1+j2+next2);

                flip_info fifi2, fifj2, fie2;
                if (check_face(c,i,fifi2))
                    fpq.push(fifi2);
                if (check_face(c,j,fifj2))
                    fpq.push(fifj2);
                if (check_edge(c,next,k,fie2))
                    fpq.push(fie2);

                nf++;
            }
        }
        else // it is an edge flip
        {
            Facet f;
            if (perform_edge_flip(fi,f))
            {
                // check the edges on the convex hull
                // except the ones on the new triangle
                
                flip_info fiel[3],fier[3];//1,fiel2,fiel3;
                CH n = f.first->neighbor(f.second);
                int in = n->index(f.first);
                for (int k = 0; k < 3; k++)
                {
                    if (check_edge(f.first,f.second,(f.second+1+k)&3,fiel[k]))
                        fpq.push(fiel[k]);
                    if (check_edge(n,in,(in+1+k)&3,fier[k]))
                        fpq.push(fier[k]);
                }
                
                
                // check faces on convex hull
                
                if (dofaces)
                {
                    flip_info fifl[3],fifr[3];//,fifl2,fifl3;
                    for (int k = 0; k < 3; k++)
                    {
                        if (check_face(f.first,(f.second+1+k)&3,fifl[k]))
                            fpq.push(fifl[k]);
                        if (check_face(n,(in+1+k)&3,fifr[k]))
                            fpq.push(fifr[k]);
                    }
                }
                ne++;
            }
        }
    }
    rt.stop();
    std::cerr << "Performed " << ne << " edge and " << nf << " face flips in " << rt.time() << " seconds" << std::endl;
    total_time += rt.time();
    std::cerr << "Total time for harmonizing Delaunay triangulation: " << total_time << " seconds" << std::endl;
}


int main(int argc, char *argv[])
{
    CLI::App app{"Harmonic triangulation"};
    
    std::string infilename = "default";
    CLI::Option *offopt = app.add_option("-o,--off", infilename, "Point set input (off)");
    CLI::Option *plyopt = app.add_option("-p,--ply", infilename, "Point set input (ply)");

    std::string delfilename = "default";
    CLI::Option *delopt = app.add_option("-d,--delaunay", delfilename, "write Delaunay tet mesh");
    std::string cgalfilename = "default";
    CLI::Option *cgalopt = app.add_option("-c,--cgaloutput", cgalfilename, "write the harmonized mesh in cgal format");
    std::string tetviewfilename = "default";
    CLI::Option *tetviewopt = app.add_option("-t,--tetviewoutput", tetviewfilename, "write the harmonized mesh in tetview format");
    bool faceflips;
    app.add_flag("-f,--faceflips", faceflips, "also check for face flips");
    
    CLI11_PARSE(app, argc, argv);

    if ( !(*offopt) && !(*plyopt))
    {
        std::cerr << "Input file required" << std::endl;
        exit(0);
    }
    std::ifstream in(infilename);
    if (!in)
    {
        std::cerr << "Error: cannot open file " << infilename << std::endl;
        exit(0);
    }
    std::vector<Point> P;
    if (*offopt)
    {
        if (!CGAL::read_off_points(in,std::back_inserter(P)))
        {
            std::cerr << "Error: cannot read off file " << infilename << std::endl;
            exit(0);
        }
    }
    if (*plyopt)
    {
        if (!CGAL::read_ply_points(in,std::back_inserter(P)))
        {
            std::cerr << "Error: cannot read ply file " << infilename << std::endl;
            exit(0);
        }
    }
    n = P.size();
    std::cout << "Read " << n << " points from file " << infilename << std::endl;

    std::vector<std::pair<Point,int> > Pi(n);
    for (int i = 0; i < n; i++)
        Pi[i] = std::make_pair(P[i],i);
    
    CGAL::Real_timer rt;
    std::cout << "Performing Delaunay triangulation on " << n << " points" << std::endl;
    rt.start();
    T = DT (Pi.begin(), Pi.end());
    rt.stop();
    std::cout << "Done in " << rt.time() << " seconds " << std::endl;
    Finite_vertices_iterator vit = T.finite_vertices_begin();
    for (int i = 0; vit != T.finite_vertices_end(); ++vit, ++i)
        vit->info() = i;
    
    if (*delopt) write_cgal(delfilename);
    stats();
    
    flip2harmonic(faceflips);
    stats();

    if (*cgalopt) write_cgal(cgalfilename);
    if (*tetviewopt) write_tetgen(tetviewfilename);

    return 0;
}
