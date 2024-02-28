#include "maintoolbar.h"
#include "ui_maintoolbar.h"
int MainToolBar::typeId = qRegisterMetaType<MainToolBar*>();
MainToolBar::MainToolBar(QWidget *parent) :
      uFunBase(parent),
      ui(new Ui::MainToolBar)
{
    ui->setupUi(this);
}

MainToolBar::~MainToolBar()
{
    delete ui;
}
