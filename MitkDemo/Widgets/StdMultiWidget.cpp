#include <QDebug>
#include <QVBoxLayout>

#include "ustatus.h"
#include "StdMultiWidget.h"
#include "ui_StdMultiWidget.h"
#include "QmitkStdMultiWidget.h"
#include "QmitkLevelWindowWidget.h"
int StdMultiWidget::typeId = qRegisterMetaType<StdMultiWidget*>();

StdMultiWidget::StdMultiWidget(QWidget *parent) :
    uFunBase(parent),
    ui(new Ui::StdMultiWidget)
{
    ui->setupUi(this);

    qDebug()<<"StdMultiWidget: create!";
}

StdMultiWidget::~StdMultiWidget()
{
    delete ui;
}

void StdMultiWidget::f_Init()
{
    qDebug()<<"StdMultiWidget: Init!";
    if (uStatus::mDataStorage == nullptr)
    {
        return ;
    }


    vlayout = new QVBoxLayout(ui->widget);
    vlayout->setMargin(0);
    vlayout->setSpacing(1);

    // Create viewParent widget with horizontal layout
    viewParent = new QWidget(ui->widget);
    vlayout->addWidget(viewParent);
    hlayout = new QHBoxLayout(viewParent);
    hlayout->setMargin(0);

    multiWidget =std::unique_ptr<QmitkStdMultiWidget>( new QmitkStdMultiWidget(viewParent));

    hlayout->addWidget(multiWidget.get());

    multiWidget->SetDataStorage(uStatus::mDataStorage);

    multiWidget->InitializeMultiWidget();

    multiWidget->AddPlanesToDataStorage();

    multiWidget->setEnabled(1);


    levelWindowWidget = std::unique_ptr<QmitkLevelWindowWidget>(new QmitkLevelWindowWidget(viewParent));
    levelWindowWidget->setStyleSheet("background:#131720;color: #80848B;font-size: 10px;");

    hlayout->addWidget(levelWindowWidget.get());

    levelWindowWidget->SetDataStorage(uStatus::mDataStorage);


    hlayout->setStretchFactor(multiWidget.get(),100);
    hlayout->setStretchFactor(levelWindowWidget.get(),1);
}

void StdMultiWidget::f_Destroy()
{
    multiWidget->RemovePlanesFromDataStorage();
    multiWidget->ActivateMenuWidget(false);

    multiWidget.release();
    multiWidget.reset();
    levelWindowWidget.release();
    levelWindowWidget.reset();
    delete hlayout;
    delete viewParent;
    delete vlayout;
}

int StdMultiWidget::f_BeforeClose()
{
    qDebug() << "StdMultiWidget::f_BeforeClose() ... ...";
    f_Destroy();
    f_Close();
    return 0;
}
void StdMultiWidget::f_Close()
{
    this->close();
}

void StdMultiWidget::f_Reset()
{
    this->multiWidget->ResetCrosshair();
}

QmitkStdMultiWidget *StdMultiWidget::GetMultiWidget()
{
    if(multiWidget)
    {
        return multiWidget.get();
    }
}
