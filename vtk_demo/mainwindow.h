#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <QHBoxLayout>
#include <QMap>
#include <itkImage.h>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkImageActor.h>
#include "LineInteractorStyle2.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
protected:
   
protected slots:
    void on_pushButton_openFile_clicked();
    void on_pushButton_openFile_2_clicked();
    void on_pushButton_show_clicked();
private:
    Ui::MainWindow *ui;


private slots:
private:
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> mrenderwindow;
    vtkSmartPointer<vtkRenderer> mrenderer;
    vtkSmartPointer<vtkImageActor> mimageactor;
    vtkSmartPointer<vtkImageActor> mimageactor2;
    vtkSmartPointer<LineInteractorStyle2> mstyle;
    vtkSmartPointer<vtkImageData> mimagedata;
    vtkSmartPointer<vtkImageData> mimagedata2;
};

#endif // MAINWINDOW_H
