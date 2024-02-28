#include "ProgressWidget.h"
#include "ui_ProgressWidget.h"

#include <QRunnable>
#include <QThread>
#include <QThreadPool>
#include <QWaitCondition>
#include <QDebug>
#include <QMovie>
#include <QKeyEvent>
#include <QPainter>
#include <QTimer>
ProgressWidget::ProgressWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ProgressWidget)
{
    qDebug()<<"ProgressWidget:current thread:"<<QThread::currentThreadId();
    ui->setupUi(this);
    ui->widget_title->hideCloseButton();
    mCanQuit = false;
}

ProgressWidget::~ProgressWidget()
{
    delete ui;
}

void ProgressWidget::Waiting(QString pTitle, QString pTips, bool pCanQuit)
{
    qDebug()<<"ProgressWidget.Waiting:current thread="<<QThread::currentThreadId();
    this->ui->widget_title->setTitle(pTitle);
    this->ui->progressBar->setVisible(false);
    this->ui->labelPercent->setVisible(false);
    this->ui->label_tips->setText(pTips);
    SetCanQuit(pCanQuit);
    this->setWindowFlags(Qt::FramelessWindowHint);
    this->show();
    QApplication::processEvents();
    qDebug()<<"ProgressWidget.Waiting finished";
}

void ProgressWidget::Start(QString pTitle, int pMax, int pMin){
    qDebug()<<"ProgressWidget.Start:current thread="<<QThread::currentThreadId();
    this->ui->progressBar->setVisible(true);
    this->ui->progressBar->setRange(pMin, pMax);
    this->ui->progressBar->setValue(pMin);
    this->ui->widget_title->setTitle(pTitle);
    this->ui->labelPercent->setVisible(true);
    this->ui->labelPercent->setText("0%");
    if(!this->windowFlags().testFlag(Qt::FramelessWindowHint))
        this->setWindowFlags(Qt::FramelessWindowHint);
    mMin = pMin;
    mMax = pMax;
    this->show();
    QApplication::processEvents();
    qDebug()<<"ProgressWidget.Start:started";
}

void ProgressWidget::SetTips(QString pTips){
    if(!pTips.isEmpty()){
        this->ui->label_tips->setText(pTips);
    }
}

void ProgressWidget::Step(int pStep, QString pTips){
//    qDebug()<<"ProgressWidget.Step:current thread="<<QThread::currentThreadId()<<"step="<<pStep;
    QString percent = QString("%1%").arg(pStep*100/mMax);
    this->ui->labelPercent->setText(percent);
    this->ui->progressBar->setValue(pStep);
    mStep = pStep;
    SetTips(pTips);
    QApplication::processEvents();
}

void ProgressWidget::Finish(){
    this->ui->labelPercent->setText("100%");
    this->ui->progressBar->setValue(mMax);
    QApplication::processEvents();
    this->hide();
}

void ProgressWidget::keyPressEvent(QKeyEvent *event)
{
    if(event->key()==Qt::Key_Escape && mCanQuit)
    {
        this->hide();
    }
}

void ProgressWidget::SetCanQuit(bool canQuit)
{
    mCanQuit = canQuit;
}
