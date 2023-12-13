#include "ReamWidget.h"
#include <QDir>
#include <QDebug>
#include <QPluginLoader>

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
#include <QDateTime>
#include <QtConcurrent/QtConcurrent>
#include <vtkSTLWriter.h>
#include <vtkWindowedSincPolyDataFilter.h>

void ReamWidget::setTool(vtkSmartPointer<vtkPolyData> source)
{
    if(!m_tool)
		m_tool = vtkSmartPointer<vtkPolyData>::New();
    m_tool->DeepCopy(source);

    mitk::DataNode::Pointer dt = m_ds->GetNamedNode("cutter");
    if(!dt)
    {
        mitk::Surface::Pointer sur = mitk::Surface::New();
        sur->SetVtkPolyData(m_tool);
        mitk::DataNode::Pointer ddt = mitk::DataNode::New();
        ddt->SetData(sur);
        ddt->SetName("cutter");
        m_ds->Add(ddt);
    }else
    {
        mitk::Surface::Pointer sur = dynamic_cast<mitk::Surface*>(dt->GetData());
        if(sur)
        {
            sur->SetVtkPolyData(m_tool);
        }
    }
    
}

void ReamWidget::setToolMatrix(vtkSmartPointer<vtkMatrix4x4> mt)
{
    if (!m_tool_matrix)
        m_tool_matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    m_tool_matrix->DeepCopy(mt);
    m_ds->GetNamedNode("cutter")->GetData()->GetGeometry()->SetIndexToWorldTransformByVtkMatrix(m_tool_matrix);
    mitk::RenderingManager::GetInstance()->RequestUpdate(m_w->GetVtkRenderWindow());
}

void ReamWidget::updateResult()
{
    if (!m_sourceImage_generated)
    {
        generateSourceImage();
        m_thread->start();
        m_sourceImage_generated = 1;
    }
    static int i = 0;
    ++i;
    Task* task = new Task(transformPolyData(m_tool_matrix, m_tool), this, i);
    QThreadPool::globalInstance()->start(task);
}

void ReamWidget::setSource(vtkSmartPointer<vtkPolyData> source)
{
    if (!m_source)
        m_source = vtkSmartPointer<vtkPolyData>::New();
    m_source->DeepCopy(source);
}

void ReamWidget::setProsthesis(vtkSmartPointer<vtkPolyData> source, vtkSmartPointer<vtkMatrix4x4> mt)
{
    if (!m_prosthesis)
        m_prosthesis = vtkSmartPointer<vtkPolyData>::New();
    if (!m_prosthesis_matrix)
        m_prosthesis_matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    m_prosthesis->DeepCopy(source);
    m_prosthesis_matrix->DeepCopy(mt);
}

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
vtkSmartPointer<vtkPolyData> ReamWidget::transformPolyData(vtkSmartPointer<vtkMatrix4x4> mt, vtkSmartPointer<vtkPolyData> p)
{
    if (!p)
        return nullptr;
    if (!mt)
        mt = vtkSmartPointer<vtkMatrix4x4>::New();
    vtkSmartPointer<vtkPolyData> ret = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkTransform> transform= vtkSmartPointer<vtkTransform>::New();
    transform->SetMatrix(mt);
    vtkSmartPointer<vtkTransformPolyDataFilter> fi= vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    fi->SetTransform(transform);
    fi->SetInputData(p);
    fi->Update();
    ret->DeepCopy(fi->GetOutput());
    return ret;
}

void ReamWidget::generateSourceImage()
{
    if (m_sourceImage_generated)
        return;
    auto prosthesis = transformPolyData(m_prosthesis_matrix, m_prosthesis);
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

    vtkSmartPointer<vtkImageStencil> stencil_femur = vtkSmartPointer<vtkImageStencil>::New();
    vtkSmartPointer<vtkPolyDataToImageStencil> pti_femur = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    vtkSmartPointer<vtkImageData> img_femur = generateImageData();
    polyDataToImageData(m_source, img_femur, pti_femur, stencil_femur);

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
    vtkSmartPointer<vtkWindowedSincPolyDataFilter> smooth = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
    smooth->SetInputData(flying->GetOutput());
    smooth->SetNumberOfIterations(20);
    smooth->Update();

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

void ReamWidget::slot_timeout()
{
    m_lock.lock();
    if(m_map.size()==0)
    {
        m_lock.unlock();
        mitk::RenderingManager::GetInstance()->RequestUpdate(m_w->GetVtkRenderWindow());
        return;
    }
    vtkSmartPointer<vtkPolyData> p= m_map.begin()->second;
    m_map.erase(m_map.begin());
    m_lock.unlock();
    mitk::Surface::Pointer sur = static_cast<mitk::Surface*>(m_ds->GetNamedNode("rest")->GetData());
    sur->SetVtkPolyData(p);
    m_ds->GetNamedNode("cutter")->GetData()->GetGeometry()->SetIndexToWorldTransformByVtkMatrix(m_tool_matrix);
    mitk::RenderingManager::GetInstance()->RequestUpdate(m_w->GetVtkRenderWindow());
}

void ReamWidget::Task::run()
{
    if (!m_polydata)
        return;
    reamwidget->m_lock.lock();
    if(reamwidget->m_map.size()>0)
    {
	    for(auto iter : reamwidget->m_map)
	    {
            if (m_order < iter.first)
            {
                reamwidget->m_lock.unlock();
                return;
            }
	    }
    }
    reamwidget->m_lock.unlock();
    vtkSmartPointer<vtkPolyDataToImageStencil> pti = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
	vtkSmartPointer<vtkImageData> img = reamwidget->generateImageData();
    vtkSmartPointer<vtkImageStencil> stencil_cutter = vtkSmartPointer<vtkImageStencil>::New();
    reamwidget->polyDataToImageData(m_polydata, img, pti, stencil_cutter);
    if (stencil_cutter->GetOutput()->GetScalarRange()[1] == 0)
        return;
    reamwidget->m_lock.lock();
#pragma omp parallel for
    for (int i = 0; i < reamwidget->mDimensions[2]; ++i) {
        for (int j = 0; j < reamwidget->mDimensions[1]; ++j) {
            for (int k = 0; k < reamwidget->mDimensions[0]; ++k) {
                uchar* pCutter = (uchar*)(stencil_cutter->GetOutput()->GetScalarPointer(k, j, i));
                uchar* pAll = (uchar*)(reamwidget->m_imageMathematicsAdd->GetOutput()->GetScalarPointer(k, j, i));
                if (*pCutter != 0) {
                    *pAll = 0;
                }
            }
        }
    }
    reamwidget->m_imageMathematicsAdd->GetOutput()->Modified();
    vtkSmartPointer<vtkDiscreteFlyingEdges3D> flying= vtkSmartPointer<vtkDiscreteFlyingEdges3D>::New();
    flying->SetInputData(reamwidget->m_imageMathematicsAdd->GetOutput());
    flying->SetNumberOfContours(3);
    flying->SetValue(0, 1);
    flying->SetValue(1, 2);
    flying->SetValue(2, 4);
    flying->SetComputeGradients(false);
    flying->SetComputeNormals(false);
    flying->SetComputeScalars(true);
    flying->Update();
    reamwidget->m_lock.unlock();
    vtkSmartPointer<vtkWindowedSincPolyDataFilter> smooth = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
    smooth->SetInputData(flying->GetOutput());
    smooth->SetNumberOfIterations(20);
    smooth->Update();
    m_polydata->DeepCopy(smooth->GetOutput());
    reamwidget->m_lock.lock();
    reamwidget->m_map.emplace(m_order,m_polydata);
    reamwidget->m_lock.unlock();
}

ReamWidget::ReamWidget(QWidget *parent)
    : WidgetBase(parent)
{
    QmitkRegisterClasses();
    vtkOutputWindow::SetGlobalWarningDisplay(0);
    QVBoxLayout* ly = new QVBoxLayout(this);
    m_w = new QmitkRenderWindow(this);
    m_w->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard3D);
    m_ds = mitk::StandaloneDataStorage::New();
    m_w->GetRenderer()->SetDataStorage(m_ds);
    ly->addWidget(m_w);
    ly->setContentsMargins(0, 0, 0, 0);

    m_thread = new QThread;
    m_timer = new QTimer;
    m_timer->setInterval(1000./30);
    m_timer->moveToThread(m_thread);
    connect(m_thread, SIGNAL(started()), m_timer, SLOT(start()));
    connect(m_timer, SIGNAL(timeout()), this, SLOT(slot_timeout()), Qt::DirectConnection);
    connect(m_thread, &QThread::finished, m_timer, &QObject::deleteLater);
    QThreadPool::globalInstance()->setMaxThreadCount(10);
}

ReamWidget::~ReamWidget()
{
    m_thread->quit();
    m_thread->wait();
}


