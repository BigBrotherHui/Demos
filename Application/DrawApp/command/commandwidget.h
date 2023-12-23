#ifndef COMMANDWIDGET_H
#define COMMANDWIDGET_H

#include <QWidget>

namespace Ui {
class CommandWidget;
}
class ReceiverBase;
class CommandWidget : public QWidget {
  Q_OBJECT
 public:
  explicit CommandWidget(QWidget *parent = nullptr);
  ~CommandWidget();
  void Register(ReceiverBase *o);
  void UnRegister(ReceiverBase *o);
  void paintEvent(QPaintEvent *event) override;

 protected:
  void notify(QString cmd);
 signals:
  void signal_sendCommand(QString cmd);
 private slots:
  void on_lineEdit_returnPressed();

 private:
  Ui::CommandWidget *ui;
  QList<ReceiverBase *> mObserver;
};

#endif  // COMMANDWIDGET_H
