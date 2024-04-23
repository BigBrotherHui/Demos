#pragma once
#include "Fitter.h"//�а���Eigenͷ�ļ���ͷ�ļ������������mitk�ļ�ǰ��
#include <QWidget>
#undef REGISTERED//����������mitkͷ�ļ�ǰ��
#include <mitkStandaloneDataStorage.h>
#include <mitkDataNode.h>
#include "PluginInterface.h"
#include <mitkDisplayActionEventBroadcast.h>
#include <mitkDisplayActionEventHandler.h>
#include "ui_widget.h"
#include <vtkPolyData.h>
class QmitkRenderWindow;
class Widget : public WidgetBase
{
    Q_OBJECT

public:
    Widget(QWidget *parent = nullptr);
    ~Widget();
    void addNode(mitk::DataNode::Pointer dt);
    void addNodeByPoint(double *pt,std::string name);
    mitk::DataNode::Pointer getNode(std::string name);
    void getPointByName(double *,std::string name);
    void addPolyData(vtkSmartPointer<vtkPolyData> p, double opa, double* color);
private slots:
    void on_pushButton_clicked();
protected:
private slots:

private:
    Ui::Widget* ui;
    QmitkRenderWindow* m_rw;
    mitk::StandaloneDataStorage::Pointer m_ds;
    mitk::DisplayActionEventBroadcast::Pointer m_DisplayActionEventBroadcast;
    std::unique_ptr<mitk::DisplayActionEventHandler> m_DisplayActionEventHandler;
};



