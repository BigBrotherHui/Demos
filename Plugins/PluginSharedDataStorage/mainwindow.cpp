#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QHBoxLayout>
#include <QmitkStdMultiWidget.h>
#include <mitkIOUtil.h>
#include <QmitkRenderWindow.h>
#include <QmitkRenderWindowWidget.h>
#include <mitkDataNode.h>
#include <mitkPointSet.h>
#include <mitkPointSetDataInteractor.h>
#include <mitkNodePredicateDataType.h>
#include <mitkImage.h>
#include <QmitkRegisterClasses.h>
//#include <mitkSplineVtkMapper3D.h>
MainWindow::MainWindow(QWidget *parent)
    : WidgetBase(parent)
    , ui(new Ui::MainWindow)
{
    QmitkRegisterClasses();
    ui->setupUi(this);
    ds1=mitk::StandaloneDataStorage::New();
    ds2=mitk::StandaloneDataStorage::New();
    ds3=mitk::StandaloneDataStorage::New();
    w1=new QmitkStdMultiWidget(this);
    w2=new QmitkStdMultiWidget(this);
    rw=new QmitkRenderWindow(this);
    QHBoxLayout *l=new QHBoxLayout(ui->widget);
    l->addWidget(w1);
    l->addWidget(w2);
    l->addWidget(rw);
    l->setStretch(0,1);
    l->setStretch(1,1);
    l->setStretch(2,1);
    w1->SetDataStorage(ds1);
    w2->SetDataStorage(ds2);
    rw->GetRenderer()->SetDataStorage(ds3);
    //rw->GetSliceNavigationController()->SetDefaultViewDirection(mitk::AnatomicalPlane::Coronal);
    rw->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard3D);
    w1->InitializeMultiWidget();
    w2->InitializeMultiWidget();
    w1->AddPlanesToDataStorage();
    w2->AddPlanesToDataStorage();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_loaddata_clicked()
{
    mitk::IOUtil::Load("D:/Images/knee-small",*ds1);
    mitk::IOUtil::Load("D:/Images/Pelvis_Right.stl", *ds2);
    auto subset=ds1->GetSubset(mitk::NodePredicateDataType::New("Image"));
    for(auto iter : *subset){
        iter->SetProperty("volumerendering",mitk::BoolProperty::New(1));
    }
    w1->ResetCrosshair();
}


void MainWindow::on_pushButton_addpoint_clicked()
{
    mitk::DataNode::Pointer dt=mitk::DataNode::New();
    mitk::PointSet::Pointer ps=mitk::PointSet::New();
    dt->SetData(ps);
    //mitk::SplineVtkMapper3D::Pointer mapper = mitk::SplineVtkMapper3D::New();
    //dt->SetMapper(mitk::BaseRenderer::Standard3D, mapper);
    dt->SetProperty("pointsize",mitk::FloatProperty::New(10.));
    mitk::PointSetDataInteractor::Pointer inter=mitk::PointSetDataInteractor::New();
    inter->SetMaxPoints(10);
    inter->LoadStateMachine("PointSet.xml");
    inter->SetEventConfig("PointSetConfig.xml");
    inter->SetDataNode(dt);
    ds1->Add(dt);
    ds2->Add(dt);
}

void MainWindow::on_pushButton_sharedTo2_clicked()
{
    w2->SetDataStorage(ds2);
    w2->ResetCrosshair();
}

void MainWindow::on_pushButton_sharedTo3_clicked()
{
    w2->SetDataStorage(ds1);
    w2->ResetCrosshair();
    /*rw->GetRenderer()->SetDataStorage(ds1);
    mitk::RenderingManager::GetInstance()->InitializeViewByBoundingObjects(rw->GetVtkRenderWindow(), ds1);*/
}
