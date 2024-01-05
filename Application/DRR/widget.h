#pragma once

#include <QWidget>
#include <QmitkRenderWindow.h>
#include <mitkStandaloneDataStorage.h>
#include <QmitkLevelWindowWidget.h>
#include <mitkDisplayActionEventBroadcast.h>
#include <mitkDisplayActionEventHandler.h>
#include <itkImage.h>
QT_BEGIN_NAMESPACE
namespace Ui { class Widget; }
QT_END_NAMESPACE

class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent=nullptr);
    ~Widget();
    
private slots:
    void on_pushButton_importImage_clicked();
    void on_horizontalSlider_translate_z_valueChanged(int v);
protected slots:
    void levelWindowChanged(const mitk::LevelWindow& levelWindow);
private:
    Ui::Widget *ui;
    QmitkRenderWindow* m_renderwindow;
    QmitkRenderWindow* m_renderwindow2d;
    mitk::StandaloneDataStorage::Pointer m_data_storage_;
    mitk::StandaloneDataStorage::Pointer m_data_storage_2d;
    QmitkLevelWindowWidget* m_lw;
    mitk::DisplayActionEventBroadcast::Pointer m_DisplayActionEventBroadcast;
    std::unique_ptr<mitk::DisplayActionEventHandler> m_DisplayActionEventHandler;
    itk::Image<short, 3>::Pointer m_image;
};

