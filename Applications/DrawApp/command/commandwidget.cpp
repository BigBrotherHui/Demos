#include "commandwidget.h"
#include <QPainter>
#include "receiverbase.h"
#include "ui_commandwidget.h"
CommandWidget::CommandWidget(QWidget *parent)
    : QWidget(parent), ui(new Ui::CommandWidget) {
  ui->setupUi(this);
}

void CommandWidget::notify(QString cmd) {
  for (int i = 0; i < mObserver.size(); i++) {
    mObserver[i]->update(cmd);
  }
}

CommandWidget::~CommandWidget() { delete ui; }

void CommandWidget::Register(ReceiverBase *o) {
  if (!o || mObserver.contains(o)) return;
  mObserver.push_back(o);
}

void CommandWidget::UnRegister(ReceiverBase *o) {
  if (!o || !mObserver.contains(o)) return;
  mObserver.removeOne(o);
}

void CommandWidget::paintEvent(QPaintEvent *event) {
  QPainter painter(this);
  painter.fillRect(rect(), QBrush("#2F3443"));
}

void CommandWidget::on_lineEdit_returnPressed() {
  QString cmd = ui->lineEdit->text();
  ui->lineEdit->clear();
  ui->plainTextEdit->appendPlainText(cmd);
  emit signal_sendCommand(cmd);
}
