#pragma once

#include <QWidget>
#include <itkImage.h>
#include <vtkPolyData.h>
#include <vtkMatrix4x4.h>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkActor.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
QT_BEGIN_NAMESPACE
namespace Ui { class Widget; }
QT_END_NAMESPACE

class Widget : public QWidget
{
    Q_OBJECT
    using PixelType= float;
public:
    Widget(QWidget *parent=nullptr);
    ~Widget();
    void requestUpdateAll();
protected:
    void setImage(bool front, itk::Image<unsigned char, 3>::Pointer img);
    vtkSmartPointer<vtkPolyData> transformPolyData(vtkSmartPointer<vtkMatrix4x4> mt, vtkSmartPointer<vtkPolyData> p);
private slots:
    void on_pushButton_importImage_clicked();
    void on_horizontalSlider_scd_valueChanged(int v);
    void on_horizontalSlider_translate_x_valueChanged(int v);
    void on_horizontalSlider_translate_y_valueChanged(int v);
    void on_horizontalSlider_translate_z_valueChanged(int v);
    void on_horizontalSlider_rotate_x_valueChanged(int v);
    void on_horizontalSlider_rotate_y_valueChanged(int v);
    void on_horizontalSlider_rotate_z_valueChanged(int v);

private:
    Ui::Widget *ui;
    QVTKOpenGLNativeWidget* m_renderwindow2dfront;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> rwfront = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkSmartPointer<vtkRenderer> rfront = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkImageActor> actorfront = vtkSmartPointer<vtkImageActor>::New();
    QVTKOpenGLNativeWidget* m_renderwindow2dside;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> rwside = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    vtkSmartPointer<vtkRenderer> rside = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkImageActor> actorside = vtkSmartPointer<vtkImageActor>::New();

    itk::Image<PixelType, 3>::Pointer m_image;
};

