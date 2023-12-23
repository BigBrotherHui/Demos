#ifndef GRAPHICSINSTANCEITEM_H
#define GRAPHICSINSTANCEITEM_H

#include "graphicsrectitem.h"
class GraphicsInstanceItem : public GraphicsRectItem {
 public:
  GraphicsInstanceItem(const QRect &rect, QGraphicsItem *parent = 0);
  QRectF boundingRect() const;
  QPainterPath shape() const;
  virtual bool loadFromXml(QXmlStreamReader *xml) { return 1; }
  virtual bool saveToXml(QXmlStreamWriter *xml) { return 1; }
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
  void setPixmap(QPixmap p) {
    pixmap = p;
    update();
  }
  QPixmap pixmap;
};

#endif  // GRAPHICSINSTANCEITEM_H
