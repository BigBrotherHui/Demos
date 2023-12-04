#include <QCoreApplication>
#include <vtkAutoInit.h>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
//#undef IGL_STATIC_LIBRARY
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <iostream>
#include <QDebug>
#include <QFileInfo>
#include "vtkOFFWriter.h"
#include <vtkNew.h>
#include <vtkSTLReader.h>
#include "vtkOFFReader.h"
#include <vtkSTLWriter.h>
#include <vtkCellData.h>
#include <vtkTriangleFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataWriter.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)
VTK_MODULE_INIT(vtkRenderingFreeType)
Eigen::MatrixXd VA,VB,VC;
Eigen::VectorXi J,I;
Eigen::MatrixXi FA,FB,FC;
igl::MeshBooleanType boolean_type(
  igl::MESH_BOOLEAN_TYPE_INTERSECT);

const char * MESH_BOOLEAN_TYPE_NAMES[] =
{
  "Union",
  "Intersect",
  "Minus",
  "XOR",
  "Resolve",
};
void update(igl::opengl::glfw::Viewer &viewer)
{
  igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,boolean_type,VC,FC,J);
  Eigen::MatrixXd C(FC.rows(),3);
  for(size_t f = 0;f<C.rows();f++)
  {
    if(J(f)<FA.rows())
    {
      C.row(f) = Eigen::RowVector3d(1,0,0);
    }else
    {
      C.row(f) = Eigen::RowVector3d(0,1,0);
    }
  }
  viewer.data().clear();
  viewer.data().set_mesh(VC,FC);
  viewer.data().set_colors(C);
  std::cout<<"A "<<MESH_BOOLEAN_TYPE_NAMES[boolean_type]<<" B."<<std::endl;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
  switch(key)
  {
    default:
      return false;
    case '.':
      boolean_type =
        static_cast<igl::MeshBooleanType>(
          (boolean_type+1)% igl::NUM_MESH_BOOLEAN_TYPES);
      break;
    case ',':
      boolean_type =
        static_cast<igl::MeshBooleanType>(
          (boolean_type+igl::NUM_MESH_BOOLEAN_TYPES-1)%
          igl::NUM_MESH_BOOLEAN_TYPES);
      break;
    case '[':
      viewer.core().camera_dnear -= 0.1;
      return true;
    case ']':
      viewer.core().camera_dnear += 0.1;
      return true;
  }
  update(viewer);
  return true;
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc,argv);
  using namespace Eigen;
  using namespace std;
  vtkNew<vtkSTLReader> reader;
  reader->SetFileName("1.stl");
  reader->Update();
  vtkNew<vtkOFFWriter> writer;
  writer->SetInputData(reader->GetOutput());
  writer->SetFileName("1.off");
  writer->Update();
  writer->Write();
  vtkNew<vtkSTLReader> reader2;
  reader2->SetFileName("2.stl");
  reader2->Update();
  vtkNew<vtkOFFWriter> writer2;
  writer2->SetInputData(reader2->GetOutput());
  writer2->SetFileName("2.off");
  writer2->Update();
  writer2->Write();
  igl::readOFF("1.off",VA,FA);
  igl::readOFF("2.off",VB,FB);
//  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
//  viewer.data().set_mesh(VC,FC);
//  // Initialize
 update(viewer);

  viewer.data().show_lines = true;
  viewer.callback_key_down = &key_down;
  viewer.core().camera_dnear = 3.9;


////  cout<<
////    "Press '.' to switch to next boolean operation type."<<endl<<
////    "Press ',' to switch to previous boolean operation type."<<endl<<
////    "Press ']' to push near cutting plane away from camera."<<endl<<
////    "Press '[' to pull near cutting plane closer to camera."<<endl<<
////    "Hint: investigate _inside_ the model to see orientation changes."<<endl;
  igl::writeOFF("3.stl",VC,FC);
  vtkNew<vtkOFFReader> offreader;
  offreader->SetFileName("3.off");
  offreader->Update();

//  vtkNew<vtkSurface> sur;
//  sur->CreateFromPolyData(offreader->GetOutput());
//  sur->GetCellData()->Initialize();
//  sur->GetPointData()->Initialize();
  ////  sur->DisplayMeshProperties();

  vtkNew<vtkSTLWriter> w;
  w->SetFileName("3.stl");
  w->SetInputData(offreader->GetOutput());
  w->Update();
  w->Write();
    viewer.launch();
  
  a.exec();
}
