#pragma execution_character_set("utf-8")
#include "systemtoolbar.h"
#include "ui_systemtoolbar.h"
#include <QDesktopWidget>
#include <QDateTime>
#include <QDebug>
#include "QMessageBox.h"
#include <QFileDialog>
#include <QScreen>
#include "runtimeutil.h"
SystemToolBar::SystemToolBar(QWidget *parent) :
      QWidget(parent),
      ui(new Ui::SystemToolBar)
{
    ui->setupUi(this);
}

SystemToolBar::~SystemToolBar()
{
    delete ui;
}

void SystemToolBar::on_pushButton_settings_clicked()
{
    
}

void SystemToolBar::on_pushButton_screenshot_clicked()
{
    QString fileName = QDateTime::currentDateTime().toString("yyyyMMddHHmmss") + ".bmp";
    QString selectedDir = QFileDialog::getSaveFileName(this, tr("保存截图"), "C:/Users/Administrator/Pictures/" + fileName, "BMP files (*.bmp)");
    if (selectedDir.isEmpty())
    {
        return;
    }
    QRect screenGeometry = QApplication::desktop()->screenGeometry();
    QPixmap pix = QApplication::screens().first()->grabWindow(
                QApplication::desktop()->winId(),
                0, 0,
                screenGeometry.width(),
                screenGeometry.height()
                );
    if (pix.isNull())
    {
        QMessageBox::warning(this,"提示","截屏失败！");
        return;
    }
    if (!pix.save(selectedDir))
    {
        QMessageBox::warning(this,"提示","截屏保存失败！");
        return;
    }
    QMessageBox::warning(this,"提示",QString("截屏成功！文件保存在：%1").arg(selectedDir));
}

void SystemToolBar::on_pushButton_exit_clicked()
{
    if(QMessageBox::No == QMessageBox::question(this,"提示","是否退出程序？", QMessageBox::Yes | QMessageBox::No))
        return;
    RuntimeUtil::SetRunning(false);
    qApp->exit(0);
}
