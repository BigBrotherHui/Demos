#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVertexGlyphFilter.h>

#include "vtkPolyDataMovingAverageFilter.h"
#include <vtkSTLReader.h>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <qapplication.h>
#include <qdatetime.h>
#include <qdebug.h>
void GenerateData(vtkPolyData*);

int main(int argc, char *argv[])
{
    QApplication ac(argc,argv);
    QVTKOpenGLNativeWidget w;
    vtkNew< vtkGenericOpenGLRenderWindow> renderWindow;
    w.SetRenderWindow(renderWindow);
  vtkSmartPointer<vtkPolyData> input =
    vtkSmartPointer<vtkPolyData>::New();

    std::string inputFilename = "D:/Demos/stl/out.stl";

    vtkSmartPointer<vtkSTLReader> reader =
      vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();
    input->ShallowCopy(reader->GetOutput());

    qDebug() <<"start smooth:"<< QDateTime::currentDateTime();
  // Compute the best fit plane.
  vtkSmartPointer<vtkPolyDataMovingAverageFilter> movingAverageFilter =
    vtkSmartPointer<vtkPolyDataMovingAverageFilter>::New();
  movingAverageFilter->SetRadius(.4);
  movingAverageFilter->SetInputData(input);
  movingAverageFilter->Update();
  qDebug() << "end smooth:" << QDateTime::currentDateTime();
  // Define viewport ranges
  // (xmin, ymin, xmax, ymax)
  double leftViewport[4] = {0.0, 0.0, 0.5, 1.0};
  double rightViewport[4] = {0.5, 0.0, 1.0, 1.0};

  // Create a mapper and actor
  vtkSmartPointer<vtkPolyDataMapper> inputMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  inputMapper->SetInputData(input);

  vtkSmartPointer<vtkActor> inputActor =
    vtkSmartPointer<vtkActor>::New();
  inputActor->SetMapper(inputMapper);

  vtkSmartPointer<vtkPolyDataMapper> averagedMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  averagedMapper->SetInputConnection(movingAverageFilter->GetOutputPort());

  vtkSmartPointer<vtkActor> averagedActor =
    vtkSmartPointer<vtkActor>::New();
  averagedActor->SetMapper(averagedMapper);

  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> leftRenderer =
    vtkSmartPointer<vtkRenderer>::New();
  leftRenderer->SetViewport(leftViewport);

  vtkSmartPointer<vtkRenderer> rightRenderer =
    vtkSmartPointer<vtkRenderer>::New();
  rightRenderer->SetViewport(rightViewport);

  renderWindow->SetSize(600,300);
  renderWindow->AddRenderer(leftRenderer);
  renderWindow->AddRenderer(rightRenderer);

  // Add the actor to the scene
  leftRenderer->AddActor(inputActor);
  rightRenderer->AddActor(averagedActor);
  //renderer->SetBackground(.3, .6, .3); // Background color green

  // Render and interact
  w.show();
  return ac.exec();
}

void GenerateData(vtkPolyData* input)
{
  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->SetThetaResolution(15);
  sphereSource->SetPhiResolution(15);
  sphereSource->Update();

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  points->ShallowCopy(sphereSource->GetOutput()->GetPoints());

  double noiseMagnitude = .1;

  for(vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
    {
    double noise[3] = {vtkMath::Random(-noiseMagnitude, noiseMagnitude),
                      vtkMath::Random(-noiseMagnitude, noiseMagnitude),
                      vtkMath::Random(-noiseMagnitude, noiseMagnitude)};
    double p[3];
    points->GetPoint(i,p);
    vtkMath::Add(p,noise,p);
    points->SetPoint(i,p);
    }

  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);

  vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputData(polydata);
  glyphFilter->Update();

  input->ShallowCopy(glyphFilter->GetOutput());

}