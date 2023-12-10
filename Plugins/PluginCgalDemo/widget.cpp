﻿#include "widget.h"
#include "ui_widget.h"
#include <QDir>
#include <QDebug>
#include <QPluginLoader>
#include "PluginInterface.h"
#include <iostream>
#include <QBoxLayout>

#include <vtkOutputWindow.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polyhedron_3.h>
#include "vtkCGALPolyDataAlgorithm.h"
#include <vtkSTLReader.h>
#include <vtkPlane.h>
#include <vtkClipPolyData.h>
#include <QmitkRenderWindow.h>
#include <QmitkRegisterClasses.h>
#include <mitkDataNode.h>
#include <mitkStandaloneDataStorage.h>
#include <mitkSurface.h>//一定要将mitk的头文件放在最下面，否则会引发eigen的错误,头文件中也不能包含mitk的文件
typedef CGAL::Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                              Point_3;
typedef CGAL::Surface_mesh<Point_3>                  Surface_mesh;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
namespace SMS = CGAL::Surface_mesh_simplification;
typedef Polyhedron::Vertex_handle                             Vertex_handle;
typedef Polyhedron::Halfedge_handle                           Halfedge_handle;
typedef Polyhedron::Facet_handle                              Facet_handle;
using Graph_halfedge = boost::graph_traits<CGAL_Surface>::halfedge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

class Widget::Impl {
public:
    QmitkRenderWindow* m_w{ nullptr };
    mitk::StandaloneDataStorage::Pointer m_ds;
};
Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget),m_impl(std::make_unique<Impl>())
{
    QmitkRegisterClasses();
    vtkOutputWindow::SetGlobalWarningDisplay(0);
    ui->setupUi(this);
    QVBoxLayout* ly = new QVBoxLayout(ui->widget);
    m_impl->m_w = new QmitkRenderWindow(this);
    m_impl->m_w->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard3D);
    m_impl->m_ds = mitk::StandaloneDataStorage::New();
    m_impl->m_w->GetRenderer()->SetDataStorage(m_impl->m_ds);
    ly->addWidget(m_impl->m_w);
    connect(ui->pushButton_simplfily, &QPushButton::clicked, this, &Widget::slot_btn_clicked);
    connect(ui->pushButton_cut, &QPushButton::clicked, this, &Widget::slot_cut_clicked);
    connect(ui->pushButton_fillhole, &QPushButton::clicked, this, &Widget::slot_fillhole_clicked);
}

Widget::~Widget()
{
    delete ui;
}

void Widget::slot_cut_clicked()
{
    vtkPolyData* vtp = static_cast<mitk::Surface*>(m_impl->m_ds->GetNamedNode("node")->GetData())->GetVtkPolyData();
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(vtp->GetCenter());
    plane->SetNormal(0, 0,1);
    vtkSmartPointer< vtkClipPolyData> clipper = vtkSmartPointer< vtkClipPolyData>::New();
    clipper->SetInputData(vtp);
    clipper->SetClipFunction(plane);
    clipper->Update();
    static_cast<mitk::Surface *>(m_impl->m_ds->GetNamedNode("node")->GetData())->SetVtkPolyData(clipper->GetOutput());
    mitk::RenderingManager::GetInstance()->RequestUpdate(m_impl->m_w->GetVtkRenderWindow());
}

void Widget::slot_fillhole_clicked()
{
    vtkPolyData* vtp = static_cast<mitk::Surface*>(m_impl->m_ds->GetNamedNode("node")->GetData())->GetVtkPolyData();
    std::unique_ptr<CGAL_Mesh> cgalMesh = vtkCGALPolyDataAlgorithm::toCGAL(vtp);
    bool success = true;

    // CGAL Processing
    // ---------------

    std::vector<Graph_Verts> patch_vertices;
    std::vector<Graph_Faces> patch_facets;

    try
    {
        // collect one halfedge per boundary cycle
        std::vector<Graph_halfedge> borderCycles;
        PMP::extract_boundary_cycles(cgalMesh->surface, std::back_inserter(borderCycles));

        // fill boundary cycles
        for (Graph_halfedge h : borderCycles)
        {
            std::get<0>(PMP::triangulate_and_refine_hole(cgalMesh->surface, h,
                std::back_inserter(patch_facets), std::back_inserter(patch_vertices),
                PMP::parameters::fairing_continuity(1)));
            //triangulate_refine_and_fair_hole以弧形补面
        }
    }
    catch (std::exception& e)
    {
        return;
    }
    auto ret=vtkCGALPolyDataAlgorithm::toVTK(cgalMesh.get());
    static_cast<mitk::Surface*>(m_impl->m_ds->GetNamedNode("node")->GetData())->SetVtkPolyData(ret);
    mitk::RenderingManager::GetInstance()->RequestUpdate(m_impl->m_w->GetVtkRenderWindow());
}

void Widget::slot_btn_clicked()
{
    std::cout << "clicked" << std::endl;
    const std::string stl_file_name = "D:/kasystem/dependence/cfg/prosthesisdata/companyA/BrandTHA/Cup_56.STL";
    vtkSmartPointer<vtkSTLReader> reader =
        vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(stl_file_name.c_str());
    reader->Update();
    std::unique_ptr<CGAL_Mesh> CGAL_Surface= vtkCGALPolyDataAlgorithm::toCGAL(reader->GetOutput());
    if (!CGAL::is_triangle_mesh(CGAL_Surface->surface))
    {
        std::cerr << "Input geometry is not triangulated." << std::endl;
        return;
    }

    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

    // In this example, the simplification stops when the number of undirected edges
    // drops below 10% of the initial count
    double stop_ratio = 0.1;
    SMS::Count_ratio_stop_predicate<Surface_mesh> stop(stop_ratio);

    int r = SMS::edge_collapse(CGAL_Surface->surface, stop);

    std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

    //std::cout << "\nFinished!\n" << r << " edges removed.\n" << surface_mesh.number_of_edges() << " final edges.\n";
    std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms" << std::endl;
    vtkSmartPointer<vtkPolyData> v=vtkCGALPolyDataAlgorithm::toVTK(CGAL_Surface.get());
    mitk::Surface::Pointer sur = mitk::Surface::New();
    sur->SetVtkPolyData(v);
    mitk::DataNode::Pointer dt = mitk::DataNode::New();
    dt->SetData(sur);
    dt->SetName("node");
    m_impl->m_ds->Add(dt);
    mitk::RenderingManager::GetInstance()->InitializeViewByBoundingObjects(m_impl->m_w->GetVtkRenderWindow(), m_impl->m_ds);
}