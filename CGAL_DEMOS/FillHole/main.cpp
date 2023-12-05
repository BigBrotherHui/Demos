#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
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
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef Surface_mesh::Vertex_index vertex_descriptor;
typedef Surface_mesh::Face_index face_descriptor;
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
int main()
{
    //构造平面
    //CGAL::Plane_3 cut_plane(0, 0, 1, 0);

    //读取模型
	std::string stl_file_name = "D:/Demos_BUILD/lib/x64/Release/stl/Femur_left.stl";
	Polyhedron_3 poly_Partition;
	CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(stl_file_name, poly_Partition);
    //PMP::clip(poly_Partition, cut_plane, CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, poly_Partition)));
    PMP::stitch_borders(poly_Partition, CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, poly_Partition)));
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
