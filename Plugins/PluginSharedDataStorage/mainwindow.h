#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <mitkStandaloneDataStorage.h>
#include "PluginInterface.h"
QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE
class QmitkStdMultiWidget;
class QmitkRenderWindow;
class MainWindow : public WidgetBase
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
protected:
    void processEvent(const PluginMetaData&) override;
private slots:
    void on_pushButton_loaddata_clicked();
    void on_pushButton_addpoint_clicked();
    void on_pushButton_sharedTo2_clicked();
    void on_pushButton_sharedTo3_clicked();
private:
    Ui::MainWindow *ui;
    QmitkStdMultiWidget *w1;
    QmitkStdMultiWidget *w2;
    mitk::StandaloneDataStorage::Pointer ds1;
    mitk::StandaloneDataStorage::Pointer ds2;
    mitk::StandaloneDataStorage::Pointer ds3;
    QmitkRenderWindow *rw;
};
#endif // MAINWINDOW_H
