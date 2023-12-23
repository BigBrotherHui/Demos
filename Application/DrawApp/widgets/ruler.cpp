#include "ruler.h"
#include <QPainter>
#include "ui_ruler.h"
Ruler::Ruler(QWidget *parent) : QWidget(parent), ui(new Ui::Ruler) {
  ui->setupUi(this);
  setAttribute(Qt::WA_TranslucentBackground);
  setAttribute(Qt::WA_TransparentForMouseEvents);
  ui->progressBar->setMinimum(0);
  ui->progressBar->setMaximum(500);
  ui->progressBar->setValue(250);
}

Ruler::~Ruler() { delete ui; }

void Ruler::updateValue(int v) {
  ui->progressBar->setValue(ui->progressBar->value() + v);
  ui->label->setText(
      QString("%1um").arg(QString::number(ui->progressBar->value() / 100.)));
  emit signal_zoomFactorChanged(ui->progressBar->value());
}
