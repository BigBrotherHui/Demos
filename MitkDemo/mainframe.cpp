#include "login.h"
#include "mainframe.h"
#include "ui_mainframe.h"
#include "ustatus.h"
#include <QDebug>

int MainFrame::typeId = qRegisterMetaType<MainFrame*>();
MainFrame::MainFrame(QWidget *parent) :
      uFunBase(parent),
      ui(new Ui::MainFrame)
{
    ui->setupUi(this);
    setWindowFlag(Qt::FramelessWindowHint);
    setWindowTitle("MitkDemo");
    if(!mLogin)
    {
        mLogin=new Login(this);
        ui->stackedWidget->addWidget(mLogin);
    }
    ui->stackedWidget->setCurrentWidget(mLogin);
    connect(mLogin,&Login::signal_loginSucceed,this,[&]
    {
        if(!uStatus::mMain)
        {
            uStatus::mMain =new MainWindow(this);
            ui->stackedWidget->addWidget(uStatus::mMain);
        }
        ui->stackedWidget->setCurrentWidget(uStatus::mMain);
    });
}

MainFrame::~MainFrame()
{
    delete ui;
}

Login *MainFrame::GetLogin()
{
    return mLogin;
}
