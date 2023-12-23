#ifndef GLOBALSIGNALINSTANCE_H
#define GLOBALSIGNALINSTANCE_H

#include <QObject>
#include "graphicsinstanceitem.h"
class GlobalSignalInstance : public QObject {
  Q_OBJECT
 public:
  GlobalSignalInstance(QObject* parent = nullptr);
  static GlobalSignalInstance* instance();
 signals:
  void signal_wrapper_instanceItemDoubleClicked(GraphicsInstanceItem*);
  void signal_wrapper_refreshItem(GraphicsItem *);
};

#endif  // GLOBALSIGNALINSTANCE_H
