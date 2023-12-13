#pragma once

#include <QWidget>
#include <vtkImageData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkNew.h>
#include <vtkSTLReader.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkImageMathematics.h>
#include <vtkDiscreteFlyingEdges3D.h>
#include <mitkStandaloneDataStorage.h>
#include <qdebug>
#include <QThreadPool>
#include <QQueue>
#include <vtkWindowedSincPolyDataFilter.h>
#include "PluginInterface.h"
#include "ReamWidget.h"
class QmitkRenderWindow;
class PluginManager;
class TimerThread;
class QThread;
class ReamWidgetTest : public WidgetBase
{
    Q_OBJECT

public:
    ReamWidgetTest(QWidget *parent = nullptr);
    ~ReamWidgetTest();
    
protected:
    
protected slots:
    void slot_btn_clicked();
private:
    ReamWidget* m_reamwidget;
    QTimer* m_timer;
};
