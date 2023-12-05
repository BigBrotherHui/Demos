#include <vtkDiscreteFlyingEdges3D.h>
#include <qapplication>
#include <vtkInformation.h>
#include <vtkImageData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkSTLReader.h>
#include <QmitkRenderWindow.h>
#include <QmitkRegisterClasses.h>
#include <mitkStandaloneDataStorage.h>
#include <mitkSurface.h>
#include <mitkDataNode.h>
#include "CustomSurfaceVtkMapper3D.h"
#include <vtkSmoothPolyDataFilter.h>
#include <vtkImageMathematics.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkOutputWindow.h>
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
vtkNew<vtkMatrix4x4> vmt_femur_prosthesis;
double mt[16]{
1.00,-0.03,-0.02,-21.90,
0.03,1.00,-0.01,-202.45,
0.02,0.01,1.00,1354.64,
0.00,0.00,0.00,1.00 };

vtkNew<vtkMatrix4x4> vmt_tibia_prosthesis;
double mt_tibia[16]{
1.00,0.02,0.00,-15.74,
- 0.02,1.00,-0.01,-208.05,
- 0.00,0.01,1.00,1325.44,
0.00,0.00,0.00,1.00
};
vtkSmartPointer<vtkPolyData> transformPolyData(vtkMatrix4x4 *mt,vtkPolyData *p)
{
    vtkSmartPointer<vtkPolyData> ret = vtkSmartPointer<vtkPolyData>::New();
    vtkNew<vtkTransform> transform;
    transform->SetMatrix(mt);
    vtkNew<vtkTransformPolyDataFilter> fi;
    fi->SetTransform(transform);
    fi->SetInputData(p);
    fi->Update();
    ret->DeepCopy(fi->GetOutput());
    return ret;
}

int main(int a, char*c[])
{
    QmitkRegisterClasses();
    vtkOutputWindow::SetGlobalWarningDisplay(0);
    QApplication app(a, c);
    QmitkRenderWindow w;
    w.GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard3D);
    mitk::StandaloneDataStorage::Pointer ds = mitk::StandaloneDataStorage::New();
    w.GetRenderer()->SetDataStorage(ds);
    
    vtkNew<vtkSTLReader> reader;
    reader->SetFileName("D:/Demos/Dependence/resources/stl/YS_NA_Tibia_PS_5_NA.STL");
    reader->Update();

    vtkNew<vtkSTLReader> reader_cutter;
    reader_cutter->SetFileName("D:/Demos/Dependence/resources/stl/thinCutter.STL");
    reader_cutter->Update();

    memcpy(vmt_tibia_prosthesis->GetData(), mt_tibia, sizeof(double) * 16);

    vtkSmartPointer<vtkImageStencil> stencil_cutter = vtkSmartPointer<vtkImageStencil>::New();
    vtkNew<vtkMatrix4x4> vmt_cutter;
    vmt_cutter->DeepCopy(vmt_tibia_prosthesis);

    auto prosthesis=transformPolyData(vmt_tibia_prosthesis, reader->GetOutput());

    double* bounds = prosthesis->GetBounds();
    mDimensions[0] = std::ceil((bounds[1] - bounds[0]) / mSpacing[0]) + 10;
    mDimensions[1] = std::ceil((bounds[3] - bounds[2]) / mSpacing[1]) + 10;
    mDimensions[2] = std::ceil((bounds[5] - bounds[4]) / mSpacing[2]) + 10;
    mExtent[1] = mDimensions[0];
    mExtent[3] = mDimensions[1];
    mExtent[5] = mDimensions[2];
    mOrigin[0] = 0.5 * mSpacing[0] + bounds[0];
    mOrigin[1] = 0.5 * mSpacing[1] + bounds[2];
    mOrigin[2] = 0.5 * mSpacing[2] + bounds[4];

    vtkSmartPointer<vtkImageStencil> stencil= vtkSmartPointer<vtkImageStencil>::New();
    vtkSmartPointer<vtkPolyDataToImageStencil> pti = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    vtkSmartPointer<vtkImageData> img = generateImageData();
    polyDataToImageData(prosthesis,img, pti, stencil);

    vtkNew<vtkSTLReader> reader_femur;
    reader_femur->SetFileName("D:/Demos/Dependence/resources/stl/TibiaLeft_Cutted.stl");
    reader_femur->Update();
    
    vtkSmartPointer<vtkImageStencil> stencil_femur = vtkSmartPointer<vtkImageStencil>::New();
    vtkSmartPointer<vtkPolyDataToImageStencil> pti_femur = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    vtkSmartPointer<vtkImageData> img_femur = generateImageData();
    polyDataToImageData(reader_femur->GetOutput(), img_femur, pti_femur, stencil_femur);

    vtkSmartPointer<vtkImageMathematics> m_imageMathematicsMultiply = vtkSmartPointer<vtkImageMathematics>::New();
    m_imageMathematicsMultiply->SetInput1Data(stencil->GetOutput());
    m_imageMathematicsMultiply->SetInput2Data(stencil_femur->GetOutput());
    m_imageMathematicsMultiply->SetOperationToMultiply();
    m_imageMathematicsMultiply->Update();

    vtkSmartPointer<vtkImageMathematics> m_imageMathematicsAdd =vtkSmartPointer<vtkImageMathematics>::New();
    m_imageMathematicsAdd->SetInput1Data(m_imageMathematicsMultiply->GetOutput());
    m_imageMathematicsAdd->SetInput2Data(stencil_femur->GetOutput());
    m_imageMathematicsAdd->SetOperationToAdd();
    m_imageMathematicsAdd->Update();

#pragma omp parallel for
    for (int i = 0; i < mDimensions[2]; ++i) {
        for (int j = 0; j < mDimensions[1]; ++j) {
            for (int k = 0; k < mDimensions[0]; ++k) {
                uchar* pMultiply = (uchar*)(m_imageMathematicsMultiply->GetOutput()->GetScalarPointer(k, j, i));
                uchar* pAll = (uchar*)(m_imageMathematicsAdd->GetOutput()->GetScalarPointer(k, j, i));
                if (*pMultiply == 1) {
                    *pAll = 2;
                }
            }
        }
    }
#pragma omp parallel for
    for (int j = 0; j < mDimensions[1]; ++j) {
        for (int k = 0; k < mDimensions[0]; ++k) {
            uchar* pAll = (uchar*)(m_imageMathematicsAdd->GetOutput()->GetScalarPointer(k, j, mDimensions[2] - 1));
            if (*pAll == 1) {
                *pAll = 4;
            }
            pAll = (uchar*)(m_imageMathematicsAdd->GetOutput()->GetScalarPointer(k, j, 0));
            if (*pAll != 0) {
                *pAll = 4;
            }
        }
    }
#pragma omp parallel for
    for (int j = 0; j < mDimensions[2]; ++j) {
        for (int k = 0; k < mDimensions[0]; ++k) {
            uchar* pAll = (uchar*)(m_imageMathematicsAdd->GetOutput()->GetScalarPointer(k, 0, j));
            if (*pAll == 1) {
                *pAll = 4;
            }
            pAll = (uchar*)(m_imageMathematicsAdd->GetOutput()->GetScalarPointer(k, mDimensions[1] - 1, j));
            if (*pAll == 1) {
                *pAll = 4;
            }
        }
    }
#pragma omp parallel for
    for (int j = 0; j < mDimensions[2]; ++j) {
        for (int k = 0; k < mDimensions[1]; ++k) {
            uchar* pAll = (uchar*)(m_imageMathematicsAdd->GetOutput()->GetScalarPointer(0, k, j));
            if (*pAll == 1) {
                *pAll = 4;
            }
            pAll = (uchar*)(m_imageMathematicsAdd->GetOutput()->GetScalarPointer(mDimensions[0] - 1, k, j));
            if (*pAll == 1) {
                *pAll = 4;
            }
        }
    }

    m_imageMathematicsAdd->GetOutput()->Modified();
    vtkNew<vtkDiscreteFlyingEdges3D> flying;
    flying->SetInputData(m_imageMathematicsAdd->GetOutput());
    flying->SetNumberOfContours(3);
    flying->SetValue(0, 1);
    flying->SetValue(1, 2);
    flying->SetValue(2, 4);
    flying->SetComputeGradients(false);
    flying->SetComputeNormals(false);
    flying->SetComputeScalars(true);
    flying->Update();

    vtkSmartPointer<vtkSmoothPolyDataFilter> smooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smooth->SetInputData(flying->GetOutput());
    smooth->SetNumberOfIterations(20);
    smooth->SetConvergence(0.0);
    smooth->SetRelaxationFactor(0.1);
    smooth->Update();

    {
        mitk::Surface::Pointer sur = mitk::Surface::New();
        sur->SetVtkPolyData(reader_cutter->GetOutput());
        mitk::DataNode::Pointer dt = mitk::DataNode::New();
        dt->SetData(sur);
        dt->SetName("cutter");
        ds->Add(dt);
        dt->GetData()->GetGeometry()->SetIndexToWorldTransformByVtkMatrix(vmt_cutter);
    }

    mitk::Surface::Pointer sur = mitk::Surface::New();
    sur->SetVtkPolyData(smooth->GetOutput());
    mitk::DataNode::Pointer dt = mitk::DataNode::New();
    dt->SetData(sur);
    ds->Add(dt);
    mitk::RenderingManager::GetInstance()->InitializeViewByBoundingObjects(w.GetVtkRenderWindow(), ds);
    dt->SetMapper(2, mitk::CustomSurfaceVtkMapper3D::New());
    w.showMaximized();

    QTimer timer;
    timer.setInterval(100);
    QObject::connect(&timer, &QTimer::timeout, [&]
        {
            static int i = 0;
            static int j = 0;
            vmt_cutter->SetElement(2, 3, 1500 - i % 400);
	        if(i>100)
		        vmt_cutter->SetElement(1, 3, vmt_cutter->GetElement(1,3)-j%20);
            i += 4;
            j += 1;
            auto r = transformPolyData(vmt_cutter, reader_cutter->GetOutput());
            vtkSmartPointer<vtkPolyDataToImageStencil> pti = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
            vtkSmartPointer<vtkImageData> img = generateImageData();
            polyDataToImageData(r, img, pti, stencil_cutter);
#pragma omp parallel for
            for (int i = 0; i < mDimensions[2]; ++i) {
                for (int j = 0; j < mDimensions[1]; ++j) {
                    for (int k = 0; k < mDimensions[0]; ++k) {
                        uchar* pCutter = (uchar*)(stencil_cutter->GetOutput()->GetScalarPointer(k, j, i));
                        uchar* pAll = (uchar*)(m_imageMathematicsAdd->GetOutput()->GetScalarPointer(k, j, i));
                        if (*pCutter != 0) {
                            *pAll = 0;
                        }
                    }
                }
            }
            m_imageMathematicsAdd->GetOutput()->Modified();
            flying->SetInputData(m_imageMathematicsAdd->GetOutput());
            flying->Update();
            smooth->SetInputData(flying->GetOutput());
            smooth->Update();
            sur->SetVtkPolyData(smooth->GetOutput());
			ds->GetNamedNode("cutter")->GetData()->GetGeometry()->SetIndexToWorldTransformByVtkMatrix(vmt_cutter);
            mitk::RenderingManager::GetInstance()->RequestUpdate(w.GetVtkRenderWindow());
        });
    timer.start();
    mitk::RenderingManager::GetInstance()->RequestUpdate(w.GetVtkRenderWindow());
    return app.exec();
}
