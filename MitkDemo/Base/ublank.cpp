#include "ublank.h"
#include "ui_ublank.h"


int uBlank::typeId = qRegisterMetaType<uBlank*>();


uBlank::uBlank(QWidget *parent) :
    uFunBase(parent),
    ui(new Ui::uBlank)
{
    ui->setupUi(this);
}

uBlank::~uBlank()
{
    delete ui;
}
