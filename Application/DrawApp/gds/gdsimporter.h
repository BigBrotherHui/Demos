#ifndef GDSIMPORTER_H
#define GDSIMPORTER_H

#include <QObject>
#include "graphicsitem.h"
class gdsimporter : public QObject {
 public:
  gdsimporter();
  static gdsimporter *instance();
  void parseFile(QString path, std::vector<GraphicsItem *> &items);
  void writeFile(QString path,const std::vector<GraphicsItem *> &items);
};

#endif  // GDSIMPORTER_H
