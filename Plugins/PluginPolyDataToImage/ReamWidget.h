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

class QmitkRenderWindow;
class PluginManager;

class ReamWidget : public QWidget
{
    Q_OBJECT

public:
    ReamWidget(QWidget *parent = nullptr);
    ~ReamWidget();
protected:
    vtkSmartPointer<vtkImageData> generateImageData();
    void polyDataToImageData(vtkSmartPointer<vtkPolyData> polydata,
        vtkSmartPointer<vtkImageData> imageData,
        vtkSmartPointer<vtkPolyDataToImageStencil> stencil,
        vtkSmartPointer<vtkImageStencil> imagestencil);
    vtkSmartPointer<vtkPolyData> transformPolyData(vtkMatrix4x4* mt, vtkPolyData* p);
protected slots:
    void slot_btn_clicked();
private slots:

private:
    double mSpacing[3]{ 0.5, 0.5, 0.5 };
    double mOrigin[3]{ 0, 0, 0 };
    int mDimensions[3]{ 160, 160, 160 };
    int mExtent[6]{ 0, 0, 0, 0, 0, 0 };
    vtkNew<vtkMatrix4x4> vmt_femur_prosthesis;
    double mt[16]{
    1.00,-0.03,-0.02,-21.90,
    0.03,1.00,-0.01,-202.45,
    0.02,0.01,1.00,1354.64,
    0.00,0.00,0.00,1.00 };

    vtkNew<vtkMatrix4x4> vmt_tibia_prosthesis;
    double mt_tibia[16]{
    1.00,0.02,0.00,-15.74,
    -0.02,1.00,-0.01,-208.05,
    -0.00,0.01,1.00,1325.44,
    0.00,0.00,0.00,1.00
    };
    vtkNew<vtkMatrix4x4> vmt_cutter;

    QmitkRenderWindow *m_w{nullptr};
    mitk::StandaloneDataStorage::Pointer m_ds;
    bool m_data_loaded{ false };
    QTimer* m_timer{nullptr};
    vtkSmartPointer<vtkSmoothPolyDataFilter> smooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    vtkNew<vtkSTLReader> reader_cutter;
    vtkSmartPointer<vtkImageStencil> stencil_cutter = vtkSmartPointer<vtkImageStencil>::New();
    vtkSmartPointer<vtkImageMathematics> m_imageMathematicsAdd = vtkSmartPointer<vtkImageMathematics>::New();
    vtkNew<vtkDiscreteFlyingEdges3D> flying;


};
#endif // WIDGET_H
