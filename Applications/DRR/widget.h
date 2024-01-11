#pragma once

#include <QWidget>
#include <QmitkRenderWindow.h>
#include <mitkStandaloneDataStorage.h>
#include <QmitkLevelWindowWidget.h>
#include <mitkDisplayActionEventBroadcast.h>
#include <mitkDisplayActionEventHandler.h>
#include <itkImage.h>
#include <vtkPolyData.h>
#include "kernel.cuh"
QT_BEGIN_NAMESPACE
namespace Ui { class Widget; }
QT_END_NAMESPACE

class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent=nullptr);
    ~Widget();
protected:
    void setImage(bool front, itk::Image<unsigned char, 3>::Pointer img);
    vtkSmartPointer<vtkPolyData> transformPolyData(vtkSmartPointer<vtkMatrix4x4> mt, vtkSmartPointer<vtkPolyData> p);
    void addPoint(double* pt);
    ;
private slots:
    void on_pushButton_importImage_clicked();
    void on_pushButton_resetView_clicked();
    void on_horizontalSlider_scd_valueChanged(int v);
    void on_horizontalSlider_translate_x_valueChanged(int v);
    void on_horizontalSlider_translate_y_valueChanged(int v);
    void on_horizontalSlider_translate_z_valueChanged(int v);
    void on_horizontalSlider_rotate_x_valueChanged(int v);
    void on_horizontalSlider_rotate_y_valueChanged(int v);
    void on_horizontalSlider_rotate_z_valueChanged(int v);
protected slots:
    void levelWindowChanged(const mitk::LevelWindow& levelWindow);
private:
    Ui::Widget *ui;
    QmitkRenderWindow* m_renderwindow;
    QmitkRenderWindow* m_renderwindow2dfront;
    QmitkRenderWindow* m_renderwindow2dside;
    mitk::StandaloneDataStorage::Pointer m_data_storage_;
    mitk::StandaloneDataStorage::Pointer m_data_storage_2dfront;
    mitk::StandaloneDataStorage::Pointer m_data_storage_2dside;
    QmitkLevelWindowWidget* m_lw;
    mitk::DisplayActionEventBroadcast::Pointer m_DisplayActionEventBroadcast;
    std::unique_ptr<mitk::DisplayActionEventHandler> m_DisplayActionEventHandler;
    itk::Image<short, 3>::Pointer m_image;
    double isocenter[3];
    vtkSmartPointer<vtkActor> m_actor_farplane;
    //vtkSmartPointer<vtkActor> actorfrustum;
    SiddonGPU* siddon{nullptr};
};

