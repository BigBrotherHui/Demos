#include "widget.h"
#include "ui_widget.h"
#include <QDir>
#include <QDebug>
#include <iostream>
#include <QmitkRenderWindow.h>
#include <mitkPointSet.h>
#include <mitkIOUtil.h>
#include <mitkDisplayActionEventHandlerStd.h>
Widget::Widget(QWidget *parent)
    : WidgetBase(parent)
{
    QVBoxLayout* l = new QVBoxLayout(this);
    m_rw = new QmitkRenderWindow(this);
    m_rw->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard2D);
    m_rw->GetRenderer()->GetSliceNavigationController()->SetDefaultViewDirection(mitk::AnatomicalPlane::Coronal);
    m_ds = mitk::StandaloneDataStorage::New();
    m_rw->GetRenderer()->SetDataStorage(m_ds);
    m_DisplayActionEventBroadcast = mitk::DisplayActionEventBroadcast::New();
    m_DisplayActionEventBroadcast->LoadStateMachine("DisplayInteraction.xml");
    m_DisplayActionEventBroadcast->SetEventConfig("DisplayConfigMITKBase.xml");
    m_DisplayActionEventBroadcast->AddEventConfig("DisplayConfigCrosshair.xml");
    m_DisplayActionEventHandler = std::make_unique<mitk::DisplayActionEventHandlerStd>();
    m_DisplayActionEventHandler->SetObservableBroadcast(m_DisplayActionEventBroadcast);
    m_DisplayActionEventHandler->InitActions();
    l->addWidget(m_rw);
    l->setContentsMargins(0, 0, 0, 0);
    //mitk::IOUtil::Load("D:/Images/knee-small", *m_ds);
    //mitk::RenderingManager::GetInstance()->InitializeViewByBoundingObjects(m_rw->GetVtkRenderWindow(), m_ds);
}

Widget::~Widget()
{
}

void Widget::addNode(mitk::DataNode::Pointer dt)
{
    if(!m_ds->Exists(dt))
    {
        m_ds->Add(dt);
        mitk::RenderingManager::GetInstance()->InitializeViewByBoundingObjects(m_rw->GetVtkRenderWindow(), m_ds);
    }
}

void Widget::addNodeByPoint(double* pt,std::string name)
{
    mitk::DataNode::Pointer dt = mitk::DataNode::New();
    mitk::PointSet::Pointer ps = mitk::PointSet::New();
    dt->SetData(ps);
    ps->SetPoint(0, mitk::Point3D(pt));
    dt->SetName(name);
    addNode(dt);
}

mitk::DataNode::Pointer Widget::getNode(std::string name)
{
    return m_ds->GetNamedNode(name);
}

void Widget::getPointByName(double *pt,std::string name)
{
    auto ps=dynamic_cast<mitk::PointSet *>(getNode(name)->GetData());
    if (!ps || ps->GetSize()==0)
       return;
    for(int i=0;i<3;++i)
    {
        pt[i] = ps->GetPoint(0)[i];
    }
}
