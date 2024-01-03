#include <vtkAutoInit.h>
#include <QApplication>
#include "ToothBite.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);
VTK_MODULE_INIT(vtkInteractionStyle);
#pragma once
#include "vtkAreaPicker.h"
#include "vtkContourFilter.h"

#include "vtkClipPolyData.h"
#include "vtkDataSetMapper.h"
#include "vtkInteractorStyleRubberBandPick.h"
#include "vtkSmartPointer.h"

#include <vtkImageCast.h>

#include <vtkInteractorStyleImage.h>
#include <vtkMetaImageReader.h>
#include <vtkSmartPointer.h>
#include "itkImageToVTKImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "vtkActor.h"
#include "vtkDecimatePro.h"
#include "vtkMarchingCubes.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkStripper.h"

#include <vtkColorTransferFunction.h>
#include <vtkGPUVolumeRayCastMapper.h>
#include <vtkObjectFactory.h>
#include <vtkPiecewiseFunction.h>
#include <vtkPlanes.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include "vtkRendererCollection.h"

#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkPolyDataToImageStencil.h>

#include <vtkAlgorithmOutput.h>
#include <vtkCutter.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkImageMapper3D.h>
#include <vtkImageStencil.h>
#include <vtkInteractorStyleImage.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkPolyDataWriter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkStripper.h>
#include <vtkVersion.h>
#include <vtkXMLPolyDataWriter.h>
#include <QDebug>

#define VTKISRBP_ORIENT 0
#define VTKISRBP_SELECT 1
class InteractorStyle
    : public vtkInteractorStyleRubberBandPick  //重载vtkInteractorStyleRubberBandPick

{
 public:
  static InteractorStyle* New();
  vtkTypeMacro(InteractorStyle, vtkInteractorStyleRubberBandPick);

  InteractorStyle() {
    selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    selectedActor = vtkSmartPointer<vtkActor>::New();
    renderer = vtkSmartPointer<vtkRenderer>::New();
  }

  virtual void OnLeftButtonUp()  //重写左键按下消息
  {
    // Forward events
    vtkInteractorStyleRubberBandPick::OnLeftButtonUp();
    if (this->CurrentMode == VTKISRBP_SELECT) {
      vtkPlanes* frustum =
          static_cast<vtkAreaPicker*>(this->GetInteractor()->GetPicker())
              ->GetFrustum();  //获得鼠标框选矩形

      vtkClipPolyData* clipper = vtkClipPolyData::New();  //裁剪polydata
      clipper->SetInputData(this->Data);
      clipper->SetClipFunction(frustum);  //!!!!很重要的一步，添加自定义隐函数
      clipper->GenerateClipScalarsOn();
      clipper->GenerateClippedOutputOn();

      clipper->SetValue(0.5);
      this->selectedMapper->SetInputConnection(clipper->GetOutputPort());
      this->selectedMapper->ScalarVisibilityOff();
      this->selectedActor->SetMapper(selectedMapper);
      this->selectedActor->GetProperty()->SetColor(1.0, 0.0, 0.0);  //(R,G,B)
      this->selectedActor->GetProperty()->SetRepresentationToWireframe();

      //	vtkSmartPointer<vtkPolyDataWriter> vtkWriter =
      // vtkSmartPointer<vtkPolyDataWriter>::New();
      // vtkWriter->SetInputConnection(clipper->GetOutputPort());
      // vtkWriter->SetFileName("D:\\test.vtk");
      // vtkWriter->Write();

      renderer->SetBackground(0.6, 0.8, 0.8);  // Blue
      renderer->AddActor(selectedActor);
      renderer->SetViewport(0.5, 0, 1, 1);

      this->Interactor->GetRenderWindow()->AddRenderer(renderer);
      //	this->GetInteractor()->GetRenderWindow()->Render();
      this->Interactor->GetRenderWindow()
          ->GetRenderers()
          ->GetFirstRenderer()
          ->Render();
    }
  }

  vtkSmartPointer<vtkRenderer> renderer;
  vtkSmartPointer<vtkPolyData> Data;
  vtkSmartPointer<vtkDataSetMapper> selectedMapper;
  vtkSmartPointer<vtkActor> selectedActor;
};
vtkStandardNewMacro(InteractorStyle);
int showVTKDataTTT(vtkSmartPointer<vtkPolyData> vtkData) {
  vtkSmartPointer<vtkContourFilter> m_skinExtractor2 =
      vtkSmartPointer<vtkContourFilter>::New();  //生成等值面/线

  m_skinExtractor2->SetInputData(vtkData);
  m_skinExtractor2->SetValue(0, 400);  //调整阈值，400以上是骨组织

  //进行精简reduce the number of triangles in a mesh
  vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
  deci->SetInputConnection(m_skinExtractor2->GetOutputPort());
  deci->SetTargetReduction(
      0.3);  // Specify the desired reduction in the total number of polygons
             // (e.g., if TargetReduction is set to 0.9, this filter will try to
             // reduce the data set to 10% of its original size).
  deci->PreserveTopologyOn();  // Turn on/off whether to preserve the topology
                               // of the original mesh. If on, mesh splitting
                               // and hole elimination will not occur. This may
                               // limit the maximum reduction that may be
                               // achieved.

  //设置优化 用拉普拉斯平滑来调整点的位置
  vtkSmartPointer<vtkSmoothPolyDataFilter> smoother =
      vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
  smoother->SetInputConnection(deci->GetOutputPort());
  smoother->SetNumberOfIterations(
      50);  // Specify the number of iterations for Laplacian smoothing

  //计算多边形网格的法线
  vtkSmartPointer<vtkPolyDataNormals> normals =
      vtkSmartPointer<vtkPolyDataNormals>::New();
  normals->SetInputConnection(smoother->GetOutputPort());
  normals->FlipNormalsOn();

  vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(normals->GetOutputPort());

  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
  actor->GetProperty()->SetColor(1.0, 1.0, 0.0);
  actor->GetProperty()->SetRepresentationToWireframe();

  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderer> renderer2 = vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  vtkSmartPointer<vtkAreaPicker> areaPicker =
      vtkSmartPointer<vtkAreaPicker>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->SetPicker(areaPicker);
  renderWindowInteractor->Initialize();

  // Set the custom stype to use for interaction.
  vtkSmartPointer<InteractorStyle> style =
      vtkSmartPointer<InteractorStyle>::New();

  style->Data = normals->GetOutput();

  renderWindowInteractor->SetInteractorStyle(style);

  renderWindow->SetSize(640, 480);

  renderer->AddActor(actor);
  renderer->ResetCamera();
  renderer->SetViewport(0, 0, 0.5, 1);
  renderer->SetBackground(0.5, 0.5, 0.5);  // Blue
  // renderer->SetBackground(.2, .3, .4);
  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
int main(int argc, char* argv[]) {
  QApplication a(argc, argv);
  //  MainWindow w;
  //  w.show();
  ToothBite t;
  t.show();
  //  vtkNew<vtkSTLReader> r;
  //  r->SetFileName("1.stl");
  //  r->Update();
  //  //  showVTKDataTTT(r->GetOutput());
  //  vtkSmartPointer<vtkContourFilter> m_skinExtractor2 =
  //      vtkSmartPointer<vtkContourFilter>::New();  //生成等值面/线

  //  m_skinExtractor2->SetInputData(r->GetOutput());
  //  m_skinExtractor2->SetValue(0, 400);  //调整阈值，400以上是骨组织

  //  //进行精简reduce the number of triangles in a mesh
  //  vtkSmartPointer<vtkDecimatePro> deci =
  //  vtkSmartPointer<vtkDecimatePro>::New();
  //  deci->SetInputConnection(r->GetOutputPort());
  //  deci->SetTargetReduction(
  //      0.3);  // Specify the desired reduction in the total number of
  //      polygons
  //             // (e.g., if TargetReduction is set to 0.9, this filter will
  //             try to
  //             // reduce the data set to 10% of its original size).
  //  deci->PreserveTopologyOn();  // Turn on/off whether to preserve the
  //  topology
  //                               // of the original mesh. If on, mesh
  //                               splitting
  //                               // and hole elimination will not occur. This
  //                               may
  //                               // limit the maximum reduction that may be
  //                               // achieved.

  //  //设置优化 用拉普拉斯平滑来调整点的位置
  //  vtkSmartPointer<vtkSmoothPolyDataFilter> smoother =
  //      vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
  //  smoother->SetInputConnection(deci->GetOutputPort());
  //  smoother->SetNumberOfIterations(
  //      50);  // Specify the number of iterations for Laplacian smoothing

  //  //计算多边形网格的法线
  //  vtkSmartPointer<vtkPolyDataNormals> normals =
  //      vtkSmartPointer<vtkPolyDataNormals>::New();
  //  normals->SetInputConnection(smoother->GetOutputPort());
  //  normals->FlipNormalsOn();

  //  vtkSmartPointer<vtkPolyDataMapper> mapper =
  //      vtkSmartPointer<vtkPolyDataMapper>::New();
  //  mapper->SetInputConnection(normals->GetOutputPort());

  //  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  //  actor->SetMapper(mapper);
  //  actor->GetProperty()->SetColor(1.0, 1.0, 0.0);
  //  actor->GetProperty()->SetRepresentationToWireframe();

  //  vtkSmartPointer<vtkRenderer> renderer =
  //  vtkSmartPointer<vtkRenderer>::New(); vtkSmartPointer<vtkRenderer>
  //  renderer2 = vtkSmartPointer<vtkRenderer>::New();

  //  vtkSmartPointer<vtkRenderWindow> renderWindow =
  //      vtkSmartPointer<vtkRenderWindow>::New();
  //  renderWindow->AddRenderer(renderer);

  //  vtkSmartPointer<vtkAreaPicker> areaPicker =
  //      vtkSmartPointer<vtkAreaPicker>::New();
  //  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
  //      vtkSmartPointer<vtkRenderWindowInteractor>::New();
  //  renderWindowInteractor->SetRenderWindow(renderWindow);
  //  renderWindowInteractor->SetPicker(areaPicker);
  //  renderWindowInteractor->Initialize();

  //  // Set the custom stype to use for interaction.
  //  vtkSmartPointer<InteractorStyle> style =
  //      vtkSmartPointer<InteractorStyle>::New();

  //  style->Data = normals->GetOutput();

  //  renderWindowInteractor->SetInteractorStyle(style);

  //  renderWindow->SetSize(640, 480);

  //  renderer->AddActor(actor);
  //  renderer->ResetCamera();
  //  renderer->SetViewport(0, 0, 0.5, 1);
  //  renderer->SetBackground(0.5, 0.5, 0.5);  // Blue
  //  // renderer->SetBackground(.2, .3, .4);
  //  renderWindow->Render();
  //  renderWindowInteractor->Start();
  return a.exec();
}
