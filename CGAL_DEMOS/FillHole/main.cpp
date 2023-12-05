﻿#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/draw_face_graph.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>
#include <QWidget>
#include <QPushButton>
using namespace std;
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> SurfaceMesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef SurfaceMesh::Vertex_index vertex_descriptor;
typedef SurfaceMesh::Face_index face_descriptor;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef Mesh::Property_map<face_descriptor, double> Face_area_map;
#define dot(u,v)   ((u).x() * (v).x() + (u).y() * (v).y() + (u).z() * (v).z())
#define norm(v)    sqrt(dot(v, v))  // norm = length of vector
#define d(u, v)     norm(u - v)       // distance = norm of difference
float pbase_Plane(K::Point_3 P, K::Plane_3 PL, K::Point_3* B)
{
    float    sb, sn, sd;
    K::Vector_3 n = PL.orthogonal_vector();
    K::Vector_3 t_p_v0 = (P - PL.point());
    sn = -dot(n, t_p_v0);
    sd = dot(n, n);
    sb = sn / sd;
    * B = P + sb * n;
    return d(P, *B);
}
struct Compute_area
{
    double operator()(const Polyhedron::Facet& f) const {
        return K::Compute_area_3()(
            f.halfedge()->vertex()->point(),
            f.halfedge()->next()->vertex()->point(),
            f.halfedge()->opposite()->vertex()->point());
    }
};
int main()
{
    //读取模型
	std::string stl_file_name = "D:/Demos_BUILD/lib/x64/Release/stl/Femur_left.stl";
    
    Polyhedron_3 poly_Partition;
    CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(stl_file_name, poly_Partition);
    poly_Partition.normalize_border();
    SurfaceMesh surface_mesh;
    CGAL::copy_face_graph(poly_Partition,surface_mesh);
    auto visited = surface_mesh.halfedges();
    std::map<SurfaceMesh::halfedge_index, bool> mp;
    for(auto h: visited)
    {
	    if(surface_mesh.is_border(h)&&!mp[h])
	    {
		    
	    }
    }
#if defined(CGAL_TEST_SUITE)
    bool cgal_test_suite = true;
#else
    bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

    if (!cgal_test_suite)
    {
        CGAL::Qt::init_ogl_context(4, 3);
        int argc = 1;
        const char* argv[2] = { "polyhedron_viewer", nullptr };
        QApplication app(argc, const_cast<char**>(argv));
        QWidget container;
        QVBoxLayout l(&container);
        CGAL::SimpleFaceGraphViewerQt
            mainwindow(app.activeWindow(), poly_Partition, "", false);
        //mainwindow.show();
        l.addWidget(&mainwindow);
        l.addWidget(new QPushButton);
        container.show();
        app.exec();
    }

	return EXIT_SUCCESS;
}
