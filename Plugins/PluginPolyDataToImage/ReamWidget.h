#ifndef REAMWIDGET_H
#define REAMWIDGET_H

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

class QmitkRenderWindow;
class PluginManager;
class TimerThread;
class QThread;
class ReamWidget : public WidgetBase
{
    Q_OBJECT

public:
    ReamWidget(QWidget *parent = nullptr);
    ~ReamWidget();
    //股骨/胫骨
    void setSource(vtkSmartPointer<vtkPolyData>source);
    //假体及规划的位置
    void setProsthesis(vtkSmartPointer<vtkPolyData>source, vtkSmartPointer<vtkMatrix4x4> mt);
    //设置打磨工具
    void setTool(vtkSmartPointer<vtkPolyData>source);
    //设置打磨工具的位置
    void setToolMatrix(vtkSmartPointer<vtkMatrix4x4> mt);
    //更新显示结果
    void updateResult();
protected:
    vtkSmartPointer<vtkImageData> generateImageData();
    void polyDataToImageData(vtkSmartPointer<vtkPolyData> polydata,
        vtkSmartPointer<vtkImageData> imageData,
        vtkSmartPointer<vtkPolyDataToImageStencil> stencil,
        vtkSmartPointer<vtkImageStencil> imagestencil);
    vtkSmartPointer<vtkPolyData> transformPolyData(vtkSmartPointer<vtkMatrix4x4> mt, vtkSmartPointer<vtkPolyData> p);
    void generateSourceImage();
protected slots:
    void slot_timeout();
private slots:

private:
    double mSpacing[3]{ 0.5, 0.5, 0.5 };
    double mOrigin[3]{ 0, 0, 0 };
    int mDimensions[3]{ 160, 160, 160 };
    int mExtent[6]{ 0, 0, 0, 0, 0, 0 };
    double mt[16]{
    1.00,-0.03,-0.02,-21.90,
    0.03,1.00,-0.01,-202.45,
    0.02,0.01,1.00,1354.64,
    0.00,0.00,0.00,1.00 };

    double mt_tibia[16]{
    1.00,0.02,0.00,-15.74,
    -0.02,1.00,-0.01,-208.05,
    -0.00,0.01,1.00,1325.44,
    0.00,0.00,0.00,1.00
    };
    QmitkRenderWindow *m_w{nullptr};
    mitk::StandaloneDataStorage::Pointer m_ds;

    vtkSmartPointer<vtkImageStencil> stencil_cutter = vtkSmartPointer<vtkImageStencil>::New();
    vtkSmartPointer<vtkImageMathematics> m_imageMathematicsAdd = vtkSmartPointer<vtkImageMathematics>::New();
    QThread* m_thread;
    QTimer* m_timer;
    std::mutex m_lock;
    std::map<int, vtkSmartPointer<vtkPolyData>> m_map;

    class Task : public QRunnable
    {
    public:
        Task(vtkSmartPointer<vtkPolyData> p,ReamWidget *w,int order)
        {
            m_polydata = p;
            reamwidget = w;
            m_order = order;
        }
		~Task()
        {
        }
        vtkSmartPointer<vtkPolyData> m_polydata{ nullptr };
        ReamWidget* reamwidget{nullptr};
        int m_order;


    protected:
        void run() override;
    };
    vtkSmartPointer<vtkPolyData> m_source{nullptr};
    vtkSmartPointer<vtkPolyData> m_tool{ nullptr };
    vtkSmartPointer<vtkPolyData> m_prosthesis{ nullptr };
    vtkSmartPointer<vtkMatrix4x4> m_prosthesis_matrix{ nullptr };
    vtkSmartPointer<vtkMatrix4x4> m_tool_matrix{ nullptr };
    bool m_sourceImage_generated{ false };
};
#endif // WIDGET_H
