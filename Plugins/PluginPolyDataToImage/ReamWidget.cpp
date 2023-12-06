#include "ReamWidget.h"
#include <QDir>
#include <QDebug>
#include <QPluginLoader>
#include "PluginInterface.h"

#include <qapplication>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <QmitkRenderWindow.h>
#include <QmitkRegisterClasses.h>
#include <mitkSurface.h>
#include <mitkDataNode.h>
#include "CustomSurfaceVtkMapper3D.h"
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkOutputWindow.h>
#include <QBoxLayout>
#include <QPushButton>
#include <QThread>
#include "TimerThread.h"
#include <vtkCamera.h>
vtkSmartPointer<vtkImageData> ReamWidget::generateImageData()
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
void ReamWidget::polyDataToImageData(vtkSmartPointer<vtkPolyData> polydata,
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
vtkSmartPointer<vtkPolyData> ReamWidget::transformPolyData(vtkMatrix4x4* mt, vtkPolyData* p)
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

void ReamWidget::hideEvent(QHideEvent* event)
{
    //return hideEvent(event);
}

void ReamWidget::showEvent(QShowEvent* event)
{
    //return QWidget::showEvent(event);
}

void ReamWidget::slot_btn_clicked()
{
    if(!m_data_loaded)
    {
        vtkNew<vtkSTLReader> reader;
        reader->SetFileName("D:/Demos/Dependence/resources/stl/YS_NA_Tibia_PS_5_NA.STL");
        reader->Update();

        reader_cutter->SetFileName("D:/Demos/Dependence/resources/stl/thinCutter.STL");
        reader_cutter->Update();

        memcpy(vmt_tibia_prosthesis->GetData(), mt_tibia, sizeof(double) * 16);

        vmt_cutter->DeepCopy(vmt_tibia_prosthesis);

        auto prosthesis = transformPolyData(vmt_tibia_prosthesis, reader->GetOutput());

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

        vtkSmartPointer<vtkImageStencil> stencil = vtkSmartPointer<vtkImageStencil>::New();
        vtkSmartPointer<vtkPolyDataToImageStencil> pti = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
        vtkSmartPointer<vtkImageData> img = generateImageData();
        polyDataToImageData(prosthesis, img, pti, stencil);

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
        flying->SetInputData(m_imageMathematicsAdd->GetOutput());
        flying->SetNumberOfContours(3);
        flying->SetValue(0, 1);
        flying->SetValue(1, 2);
        flying->SetValue(2, 4);
        flying->SetComputeGradients(false);
        flying->SetComputeNormals(false);
        flying->SetComputeScalars(true);
        flying->Update();

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
            m_ds->Add(dt);
            dt->GetData()->GetGeometry()->SetIndexToWorldTransformByVtkMatrix(vmt_cutter);
        }

        mitk::Surface::Pointer sur = mitk::Surface::New();
        sur->SetVtkPolyData(smooth->GetOutput());
        mitk::DataNode::Pointer dt = mitk::DataNode::New();
        dt->SetData(sur);
        m_ds->Add(dt);
        dt->SetName("rest");
        mitk::RenderingManager::GetInstance()->InitializeViewByBoundingObjects(m_w->GetVtkRenderWindow(), m_ds);
        dt->SetMapper(2, mitk::CustomSurfaceVtkMapper3D::New());
        vtkCamera *camera=m_w->GetRenderer()->GetVtkRenderer()->GetActiveCamera();
        camera->SetViewUp(0, 1, 0);
        camera->SetFocalPoint(smooth->GetOutput()->GetCenter());
        camera->SetPosition(smooth->GetOutput()->GetCenter()[0], smooth->GetOutput()->GetCenter()[1], smooth->GetOutput()->GetCenter()[2] + 400);
        m_w->GetRenderer()->GetVtkRenderer()->ResetCameraClippingRange();
    }
    m_thread->start();
}

void ReamWidget::slot_timeout()
{
    static int i = 0;
    static int j = 0;
    vmt_cutter->SetElement(2, 3, 1400 - i % 50);
    if (i > 100)
        vmt_cutter->SetElement(1, 3, vmt_cutter->GetElement(1, 3) - j % 20);
    i += 1;
    j += 2;
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
    mitk::Surface::Pointer sur = static_cast<mitk::Surface*>(m_ds->GetNamedNode("rest")->GetData());
    sur->SetVtkPolyData(smooth->GetOutput());
    m_ds->GetNamedNode("cutter")->GetData()->GetGeometry()->SetIndexToWorldTransformByVtkMatrix(vmt_cutter);
    mitk::RenderingManager::GetInstance()->RequestUpdate(m_w->GetVtkRenderWindow());
}

ReamWidget::ReamWidget(QWidget *parent)
    : QWidget(parent)
{
    QmitkRegisterClasses();
    vtkOutputWindow::SetGlobalWarningDisplay(0);
    QVBoxLayout* ly = new QVBoxLayout(this);
    m_w = new QmitkRenderWindow(this);
    m_w->setAttribute(Qt::WA_TransparentForMouseEvents);
    m_w->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard3D);
    m_ds = mitk::StandaloneDataStorage::New();
    m_w->GetRenderer()->SetDataStorage(m_ds);
    ly->addWidget(m_w);
    QPushButton* btn = new QPushButton(this);
    ly->addWidget(btn);
    connect(btn, &QPushButton::clicked, this, &ReamWidget::slot_btn_clicked);
    m_thread = new QThread();
    m_timer = new QTimer;
    m_timer->setInterval(100);
    m_timer->moveToThread(m_thread);
    connect(m_thread, SIGNAL(started()), m_timer, SLOT(start()));//线程打开同时启动线程中的定时器
    connect(m_timer, SIGNAL(timeout()), this, SLOT(slot_timeout()), Qt::DirectConnection);
}

ReamWidget::~ReamWidget()
{
    m_thread->quit();
    m_thread->wait();
}