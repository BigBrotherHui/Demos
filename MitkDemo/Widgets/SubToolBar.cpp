#include "SubToolBar.h"
#include "ui_SubToolBar.h"
#include <QDebug>
#include "ufunbase.h"
#include "ufunction.h"
#include "ustatus.h"
#include "umainfunbase.h"
#include "QMessageBox"
#include "global.h"

int SubToolBar::typeId = qRegisterMetaType<SubToolBar*>();

SubToolBar::SubToolBar(QWidget *parent) :
    uFunBase(parent),
    ui(new Ui::SubToolBar)
{
    ui->setupUi(this);
}

SubToolBar::~SubToolBar()
{
    delete ui;
}

void SubToolBar::f_Refresh()
{
}

void SubToolBar::f_Init()
{
}

void SubToolBar::on_pushButton_open_clicked()
{
    uMainFunBase* mMainFunBase = uFunction::getInStance()->f_GetMain();
    mMainFunBase->f_Open_Center("CaseManage");
}
