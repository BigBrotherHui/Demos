#include "Viewer.h"
#include <QMouseEvent>
#include <QGLFunctions>
#include <CGAL/Qt/CreateOpenGLContext.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
typedef CGAL::Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                              Point_3;
typedef CGAL::Surface_mesh<Point_3>                  Surface_mesh;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
namespace SMS = CGAL::Surface_mesh_simplification;

Viewer::Viewer(QWidget* parent)
  : CGAL::QGLViewer(parent),
    m_custom_mouse(false)
{
    m_pScene = new Scene();
    setScene(m_pScene);
    setManipulatedFrame(m_pScene->manipulatedFrame());
    const std::string stl_file_name = "D:/kasystem/dependence/cfg/prosthesisdata/companyA/BrandTHA/Cup_56.STL";
    Polyhedron poly_Partition;
    CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(stl_file_name, poly_Partition);
    //Surface_mesh surface_mesh;
    //CGAL::copy_face_graph(poly_Partition, surface_mesh);

    if (!CGAL::is_triangle_mesh(poly_Partition))
    {
        std::cerr << "Input geometry is not triangulated." << std::endl;
        return;
    }

    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

    // In this example, the simplification stops when the number of undirected edges
    // drops below 10% of the initial count
    double stop_ratio = 0.01;
    SMS::Count_ratio_stop_predicate<Surface_mesh> stop(stop_ratio);

    int r = SMS::edge_collapse(poly_Partition, stop);

    std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

    //std::cout << "\nFinished!\n" << r << " edges removed.\n" << surface_mesh.number_of_edges() << " final edges.\n";
    std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms" << std::endl;
    //CGAL::copy_face_graph(surface_mesh, poly_Partition);
    m_pScene->setSurfaceMesh(&poly_Partition);
    m_pScene->update_bbox();
    const Scene::Bbox bbox = m_pScene->bbox();
    const double xmin = bbox.xmin();
    const double ymin = bbox.ymin();
    const double zmin = bbox.zmin();
    const double xmax = bbox.xmax();
    const double ymax = bbox.ymax();
    const double zmax = bbox.zmax();
    CGAL::qglviewer::Vec
        vec_min(xmin, ymin, zmin),
        vec_max(xmax, ymax, zmax);
    setSceneBoundingBox(vec_min, vec_max);
    camera()->showEntireScene();
    update();
}

void Viewer::setScene(Scene* pScene)
{
    this->m_pScene = pScene;
}

void Viewer::draw()
{
  CGAL::QGLViewer::draw();
  if(m_pScene != nullptr)
  {
      m_pScene->draw(this);
  }
}

void Viewer::initializeGL()
{
  CGAL::QGLViewer::initializeGL();
  setBackgroundColor(::Qt::white);
  //m_pScene->initGL(this);
}

void Viewer::mousePressEvent(QMouseEvent* e)
{
  if ( e->modifiers() == Qt::ControlModifier )
  {
    m_custom_mouse = true;
  }

  CGAL::QGLViewer::mousePressEvent(e);
}

void Viewer::mouseReleaseEvent(QMouseEvent* e)
{
  if ( m_custom_mouse )
  {
  }

  CGAL::QGLViewer::mouseReleaseEvent(e);
}

