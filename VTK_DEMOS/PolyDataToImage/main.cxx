#include <vtkDiscreteFlyingEdges3D.h>

#include <vtkActor.h>
#include <vtkRenderer.h>

#include <QVTKOpenGLNativeWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPolyDataMapper.h>
#include <qapplication>
#include <vtkInformation.h>
#include <vtkImageData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkSTLReader.h>

double mSpacing[3]{ 0.5, 0.5, 0.5 };
double mOrigin[3]{ 0, 0, 0 };
int mDimensions[3]{ 160, 160, 160 };
int mExtent[6]{ 0, 0, 0, 0, 0, 0 };

vtkSmartPointer<vtkImageData> generateImageData()
{
    vtkSmartPointer<vtkImageData> ret = vtkSmartPointer<vtkImageData>::New();
    vtkSmartPointer<vtkInformation> in = vtkSmartPointer<vtkInformation>::New();
    ret->SetScalarType(VTK_UNSIGNED_CHAR, in);
    ret->SetSpacing(mSpacing[0], mSpacing[1], mSpacing[2]);
    ret->SetOrigin(mOrigin[0], mOrigin[1], mOrigin[2]);
    ret->SetExtent(mExtent);
    ret->SetDimensions(mDimensions);
    ret->SetNumberOfScalarComponents(1, in);
    ret->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
    return ret;
}

void polyDataToImageData(vtkSmartPointer<vtkPolyData> polydata,
    vtkSmartPointer<vtkImageData> imageData,
    vtkSmartPointer<vtkPolyDataToImageStencil> stencil,
    vtkSmartPointer<vtkImageStencil> imagestencil)
{
    vtkIdType count = imageData->GetNumberOfPoints();
#pragma omp parallel for
    for (vtkIdType i = 0; i < count; ++i) {
        imageData->GetPointData()->GetScalars()->SetTuple1(i, 1);
    }
    unsigned char outval = 0;
    if (polydata) stencil->SetInputData(polydata);
    stencil->SetOutputOrigin(mOrigin);
    stencil->SetOutputSpacing(mSpacing);
    stencil->SetOutputWholeExtent(mExtent);
    stencil->Update();

    imagestencil->SetInputData(imageData);
    if (polydata)
        imagestencil->SetStencilData(stencil->GetOutput());
    else
        imagestencil->SetStencilConnection(stencil->GetOutputPort());
    imagestencil->ReverseStencilOff();
    imagestencil->SetBackgroundValue(outval);
    imagestencil->Update();
}

int main(int a, char*c[])
{
    QApplication app(a, c);
    QVTKOpenGLNativeWidget w;
    vtkNew<vtkGenericOpenGLRenderWindow> renWin;
    w.SetRenderWindow(renWin);
    vtkNew<vtkRenderer> aRenderer;
    renWin->AddRenderer(aRenderer);

    vtkNew<vtkSTLReader> reader;
    reader->SetFileName("D:/Demos/Dependence/resources/stl/YS_NA_Femoral_PS_5_L.STL");
    reader->Update();

    vtkSmartPointer<vtkImageStencil> stencil= vtkSmartPointer<vtkImageStencil>::New();
    vtkSmartPointer<vtkPolyDataToImageStencil> pti = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    vtkSmartPointer<vtkImageData> img = generateImageData();
    polyDataToImageData(reader->GetOutput(),img, pti, stencil);

    vtkNew<vtkActor> lens;
    vtkNew<vtkPolyDataMapper> mapper;
    lens->SetMapper(mapper);
    aRenderer->AddActor(lens);

    renWin->SetSize(640, 480);
    renWin->SetWindowName("TissueLens");
    w.show();
    return app.exec();
}
