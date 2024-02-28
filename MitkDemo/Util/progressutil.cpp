#pragma execution_character_set("utf-8")
#include "progressutil.h"
#include <QApplication>
#include <QDebug>
#include <QThread>
#include <QTime>

ProgressUtil::ProgressUtil(QObject* parent)
    :QObject(parent)
{
    //TODO 打开进度窗口时，不传parent会导致进度窗口在windows任务栏中单独存在，传了parent会出现无法屏蔽操作的情况
    qDebug()<<"ProgressUtil:current thread="<<QThread::currentThreadId();
    mProgressWidget = new ProgressWidget();
    mProgressWidget->setWindowFlag(Qt::SubWindow);
    mFinish = false;
}

ProgressUtil::~ProgressUtil(){
    if(!mFinish){
        mFinish = true;
    }
    mProgressWidget->deleteLater();
}

void ProgressUtil::Waiting(QString pTitle, QString pTips, bool pCanQuit)
{
    mProgressWidget->Waiting(pTitle, pTips, pCanQuit);
}

void ProgressUtil::Start(QString pTitle, int pMax, int pMin){
    qDebug()<<"ProgressUtil.Start"<<pTitle<<pMax<<pMin;
    mProgressWidget->Start(pTitle, pMax, pMin);
    mMax = pMax;
    mMin = pMin;
    mStep = pMin*100;
}

void ProgressUtil::SetTips(QString pTips){
    mProgressWidget->SetTips(pTips);
}

void ProgressUtil::Step(int pStep, QString pTips){
    qDebug()<<"ProgressUtil.Step:"<<pStep<<pTips;
    mStep = pStep*100;
    mProgressWidget->Step(pStep, pTips);
}

void ProgressUtil::AutoStep(std::function<void ()> pCallback, QString pTips)
{
    AutoStep(pCallback, -1, pTips);
}

void ProgressUtil::AutoStep(std::function<void ()> pCallback, int pStep, QString pTips)
{
    AutoStep(pCallback, pStep, 1, 100, pTips);
}

void ProgressUtil::AutoStep(std::function<void ()> pCallback, int pStep, int pSeconds, int pInterval, QString pTips)
{
    std::function<void* ()> callback = [=]{
        pCallback();
        return new QString("");
    };
    AutoStepWithReturn(callback, pStep, pSeconds, pInterval, pTips);
}

void *ProgressUtil::AutoStepWithReturn(std::function<void *()> pCallback, QString pTips)
{
    return AutoStepWithReturn(pCallback, -1, pTips);
}

void *ProgressUtil::AutoStepWithReturn(std::function<void *()> pCallback, int pStep, QString pTips)
{
    return AutoStepWithReturn(pCallback, pStep, 1, 100, pTips);
}

void* ProgressUtil::AutoStepWithReturn(std::function<void*()> pCallback, int pStep, int pSeconds, int pInterval, QString pTips){
    qDebug()<<"ProgressUtil.AutoStepWithReturn:"<<pStep<<pSeconds<<pInterval<<pTips;
    if(pStep<0){
        pStep = mMax;
    }
    if(pSeconds<0)
    {
        pSeconds = 1;
    }
    if(pInterval<0)
    {
        pInterval = 100;
    }
//    QThread* myThread = new QThread();
//    QThread* thread = new QThread();
//    myThread->moveToThread(thread);
//    thread->start();
//    mAutoStepFinish = false;
//    connect(this, &ProgressUtil::signal_noname, myThread, [=]() mutable{
//        qDebug()<<"ProgressUtil.AutoStep:signal_noname";
//        while(!mAutoStepFinish && !mFinish){
//            QThread::msleep(pInterval);
//            mStep = mStep + (pStep*100 - mStep) / (pSeconds*1000/pInterval/2);
//            qDebug()<<"ProgressUtil.AutoStep:"<<mStep;
//            mProgressWidget->Step(mStep/100, pTips);
// //            emit signal_step(mStep/100);
//        }
//    });
//    connect(this, &ProgressUtil::signal_step, this, [=](int pStep){
//        qDebug()<<"ProgressUtil.AutoStep:step="<<pStep;
//        mProgressWidget->Step(pStep, pTips);
//    });
//    emit signal_noname();
//    mAutoStepRet = callback();
//    QApplication::processEvents();
//    AutoStepFinish();
//    thread->quit();
//    thread->wait();
//    thread->deleteLater();
//    myThread->quit();
//    myThread->wait();
//    myThread->deleteLater();
//    return mAutoStepRet;


    QThread* myThread = new QThread();
    QThread* thread = new QThread();
    myThread->moveToThread(thread);
    thread->start();
    connect(this, &ProgressUtil::signal_noname, myThread, [=]() mutable{
        emit signal_autoStepFinish(pCallback());
    });
    connect(this, &ProgressUtil::signal_autoStepFinish, this, [=](void* ret){
        mAutoStepRet = ret;
        AutoStepFinish();
    });
    emit signal_noname();
    mAutoStepFinish = false;
    while(!mAutoStepFinish && !mFinish){
        QTime dieTime = QTime::currentTime().addMSecs(pInterval);
        while( QTime::currentTime() < dieTime )
        {
            QApplication::processEvents(QEventLoop::AllEvents, pInterval);
        }
        mStep = mStep + (pStep*100 - mStep) / (pSeconds*1000/pInterval/2);
        qDebug()<<"ProgressUtil.AutoStep:"<<mStep;
        mProgressWidget->Step(mStep/100, pTips);
    }
    mStep = pStep*100;
    mProgressWidget->Step(pStep, pTips);
    thread->quit();
    thread->wait();
    thread->deleteLater();
    myThread->quit();
    myThread->wait();
    myThread->deleteLater();
    return mAutoStepRet;
}

void ProgressUtil::AutoStepWithImage(std::function<void ()> pCallback, QString pTips)
{
    AutoStepWithImage(pCallback, -1, pTips);
}

void ProgressUtil::AutoStepWithImage(std::function<void ()> pCallback, int pStep, QString pTips)
{
    AutoStepWithImage(pCallback, pStep, 1, 100, pTips);
}

void ProgressUtil::AutoStepWithImage(std::function<void ()> pCallback, int pStep, int pSeconds, int pInterval, QString pTips)
{
    qDebug()<<"ProgressUtil.AutoStepWithImage:"<<pStep<<pSeconds<<pInterval<<pTips;
    if(pStep<0){
        pStep = mMax;
    }
    if(pSeconds<0)
    {
        pSeconds = 1;
    }
    if(pInterval<0)
    {
        pInterval = 100;
    }
    //由于多线程异步执行处理图像时，有时候会出现系统崩溃的情况（原因不明），因此，暂时不作多线程异步来实现进度处理
    pCallback();
    mProgressWidget->Step(pStep, pTips);
}

void ProgressUtil::AutoStepFinish(){
    qDebug()<<"ProgressUtil.AutoStepFinish";
    mAutoStepFinish = true;
}

void ProgressUtil::Finish(){
    qDebug()<<"ProgressUtil.Finish";
    mFinish = true;
    mProgressWidget->hide();
}

bool ProgressUtil::IsQuit()
{
    return mProgressWidget->isHidden();
}

