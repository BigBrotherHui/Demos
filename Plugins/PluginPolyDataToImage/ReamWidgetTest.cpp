#include "ReamWidgetTest.h"
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
#include <vtkCamera.h>
#include <QDateTime>
#include <QtConcurrent/QtConcurrent>
#include <vtkSTLWriter.h>
#include <vtkWindowedSincPolyDataFilter.h>


void ReamWidgetTest::slot_btn_clicked()
{
    vtkSmartPointer<vtkSTLReader> reader_femur = vtkSmartPointer<vtkSTLReader>::New();
    reader_femur->SetFileName(QString(qApp->applicationDirPath() + "/stl/FemurLeft_Cutted.stl").toStdString().c_str());
    reader_femur->Update();
    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(QString(qApp->applicationDirPath() + "/stl/YS_NA_Femoral_PS_5_L.STL").toStdString().c_str());
    reader->Update();
    vtkNew<vtkSTLReader> reader_cutter;
    reader_cutter->SetFileName(QString(qApp->applicationDirPath() + "/stl/roughCutter.STL").toStdString().c_str());
    reader_cutter->Update();
    double mt_tibia[16]{
		1.00,0.02,0.00,-28.036,
		-0.02,1.00,-0.01,-204.929,
		-0.00,0.01,1.00,1355.258,
		0.00,0.00,0.00,1.00
    };
    vtkSmartPointer<vtkMatrix4x4> prosthesis_matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    memcpy(prosthesis_matrix->GetData(), mt_tibia, sizeof(double) * 16);
    m_reamwidget->setSource(reader_femur->GetOutput());
    m_reamwidget->setTool(reader_cutter->GetOutput());
    m_reamwidget->setProsthesis(reader->GetOutput(), prosthesis_matrix);
    m_reamwidget->updateResult();
    m_timer = new QTimer(this);
    m_timer->setInterval(1000 / 30);
    connect(m_timer, &QTimer::timeout, this, [&]
        {
            vtkSmartPointer<vtkMatrix4x4> vmt_cutter = vtkSmartPointer<vtkMatrix4x4>::New();
            Eigen::Vector3d v1{ -15.0517, -216.174, 1322.04 }, v2{ -24.1781, -210.9, 1323.7 };
			Eigen::Vector3d dir = (v2 - v1).normalized();
            static int ii = 0;
            ++ii;

            double pt[3]{ -15.0517, -216.174, 1322.04 };
            for (int i = 0; i < 3; ++i)
                vmt_cutter->SetElement(i, 3, pt[i]);
        
			for (int i = 0; i < 3; ++i)
			{
			    vmt_cutter->SetElement(i, 3, vmt_cutter->GetElement(i, 3) + .1*ii*dir[i]);
			}
			
            m_reamwidget->setToolMatrix(vmt_cutter);
            m_reamwidget->updateResult();
        });
    m_timer->start();
}

ReamWidgetTest::ReamWidgetTest(QWidget *parent)
    : WidgetBase(parent)
{
    QVBoxLayout* ly = new QVBoxLayout(this);
    m_reamwidget = new ReamWidget(this);
    ly->addWidget(m_reamwidget);
    QPushButton* btn = new QPushButton(this);
    connect(btn, &QPushButton::clicked, this, &ReamWidgetTest::slot_btn_clicked);
    ly->addWidget(btn);
    ly->setContentsMargins(0, 0, 0, 0);
}

ReamWidgetTest::~ReamWidgetTest()
{
}