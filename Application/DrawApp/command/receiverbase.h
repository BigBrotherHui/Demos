#ifndef RECEIVERBASE_H
#define RECEIVERBASE_H

#include <QObject>
class ReceiverBase : public QObject {
  Q_OBJECT
 public:
  ReceiverBase(QObject *parent = nullptr);
  virtual void update(QString cmd) = 0;
};

#endif  // RECEIVERBASE_H
