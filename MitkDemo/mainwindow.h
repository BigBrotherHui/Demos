
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "umainfunbase.h"
#include <QtCore/qobjectdefs.h>
#include <QtWidgets/QMainWindow>

class QWidget;
class QCloseEvent;
QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow,uMainFunBase
{
    Q_OBJECT

public:
    static int typeId;
    Q_INVOKABLE explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void f_Init();
    void setPatientInfo(QString t);
protected:
    virtual void closeEvent(QCloseEvent* e) override;
private slots:
    void slot_timeout();
private:
    Ui::MainWindow *ui;
    QTimer* mTimer;

};
#endif // MAINWINDOW_H
