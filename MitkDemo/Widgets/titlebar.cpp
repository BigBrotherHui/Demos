#include "titlebar.h"
#include "ui_titlebar.h"
#include <QPainter>
#include <QStyleOption>
#include <QDebug>
#include <QDynamicPropertyChangeEvent>
#include <QDesktopWidget>
TitleBar::TitleBar(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::TitleBar)
{
    ui->setupUi(this);
    if(!target)
        target=parent;
    take_over();
    setTarget(target);
    ui->label_title->setStyleSheet("font-size: 24px;font-family: \"思源黑体 CN Bold\";color: #FFFFFF;");
    ui->pushButton_closeWindow->setStyleSheet("border:none;background:transparent;");
}

TitleBar::~TitleBar()
{
    delete ui;
}

void TitleBar::setTitle(const QString &title)
{
    ui->label_title->setText(title);
}

QString TitleBar::getTitle()
{
    return ui->label_title->text();
}

void TitleBar::paintEvent(QPaintEvent *event)
{
    QStyleOption opt;
    opt.init(this);
    QPainter painter(this);
    style()->drawPrimitive(QStyle::PE_Widget,&opt,&painter,this);
}

bool TitleBar::event(QEvent *event)
{
    if(event->type()==QEvent::DynamicPropertyChange)
    {
        QDynamicPropertyChangeEvent *ev=static_cast<QDynamicPropertyChangeEvent *>(event);
        if(ev->propertyName()!="margin-left")
            return false;
        int marginLeft=property(ev->propertyName()).toInt();
        emit signal_propertyChanged(marginLeft);
    }
    else if(event->type()==QEvent::Show)
    {
        QDesktopWidget *deskwidget=qApp->desktop();
        int left=deskwidget->width()-target->width();
        int top=deskwidget->height()-target->height();
        left/=2;
        top/=2;
        target->move(left,top);
    }
    return QWidget::event(event);
}

void TitleBar::setTarget(QWidget *target)
{
    if(!target)
        target=static_cast<QWidget *>(parent());
    this->target=target;
    emit signal_targetChanged(target);
}

void TitleBar::hideCloseButton()
{
    ui->pushButton_closeWindow->hide();
}

void TitleBar::closeByOthers()
{
    disconnect(ui->pushButton_closeWindow,&QPushButton::clicked,0,0);
    connect(ui->pushButton_closeWindow,&QPushButton::clicked,this,&TitleBar::signal_closeWindow);
}

void TitleBar::take_over()
{
    connect(this,&TitleBar::signal_propertyChanged,this,&TitleBar::slot_propertyChanged);
    connect(this,&TitleBar::signal_targetChanged,this,&TitleBar::slot_targetChanged);
}

void TitleBar::slot_targetChanged(QWidget *target)
{
    if(!target)
        return;
    setTitle(target->windowTitle());
    target->setWindowFlags(target->windowFlags() | Qt::FramelessWindowHint);
    disconnect(ui->pushButton_closeWindow,&QPushButton::clicked,0,0);
    connect(ui->pushButton_closeWindow,&QPushButton::clicked,target,&QWidget::close);
    connect(target,&QWidget::windowTitleChanged,this,&TitleBar::setTitle);
}

void TitleBar::slot_propertyChanged(int value)
{
    if( ui && ui->horizontalLayout)
    {
        QMargins margins=ui->horizontalLayout->contentsMargins();
        ui->horizontalLayout->setContentsMargins(value,margins.top(),margins.right(),margins.bottom());
    }
}
