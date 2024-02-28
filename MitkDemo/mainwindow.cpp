#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include "ufunbase.h"
#include "ufunction.h"
#include "ustatus.h"
#include <QPointer>
#include "QMessageBox"
#include <QCloseEvent>
#include "runtimeutil.h"
#include "systemtoolbar.h"
#include <QGraphicsBlurEffect>
#include "global.h"
#include "SqlDemo.h"
#include <QTimer>
int MainWindow::typeId = qRegisterMetaType<MainWindow*>();

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //ui->widget->setVisible(false);
    f_Init();
    RuntimeUtil::SetRunning(true);
    mTimer = new QTimer(this);
    connect(mTimer,&QTimer::timeout,this,&MainWindow::slot_timeout);
    mTimer->start();
}

void MainWindow::slot_timeout()
{
    QString time = QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss");
    ui->label_time->setText(time);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::f_Init()
{
    mStackedWidget_Center = this->ui->stackedWidget;
    mStackedWidget_control_right = this->ui->stackedWidget_control_right;
    mStackedWidget_control_left = this->ui->stackedWidget_control_left;
    mStackedWidget_control_toolbar = this->ui->stackedWidget_topbar;
    mStackedWidget_control_footbar = this->ui->stackedWidget_foot;

    uFunction::getInStance()->f_InsertMain(this);

    f_Open_Center (uStatus::mMainCenterWidget);
    f_OpenControl_right(uStatus::mMainControlWidgetRight);
    f_OpenControl_toobar(uStatus::mMainToolbarWidget);
    f_OpenControl_left(uStatus::mMainControlWidgetLeft);
    f_OpenControl_footbar(uStatus::mMainFootbarWidget);
    
}

void MainWindow::setPatientInfo(QString t)
{
    ui->label_patientInfo->setText(t);
}

void MainWindow::closeEvent(QCloseEvent *e)
{
    if (QMessageBox::No == QMessageBox::question(this, tr("提示"), tr("请确认是否退出？"), QMessageBox::Yes | QMessageBox::No))
    {
        e->ignore();
        return;
    }
    e->accept();
    RuntimeUtil::SetRunning(false);
}




