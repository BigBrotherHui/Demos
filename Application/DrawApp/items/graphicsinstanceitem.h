#ifndef GRAPHICSINSTANCEITEM_H
#define GRAPHICSINSTANCEITEM_H

#include "graphicsrectitem.h"
class DrawView;
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
  QGraphicsItem *duplicate() const;
  QPixmap pixmap;
  DrawView *view{nullptr};
//  QList<GraphicsInstanceItem *> instanceItems;
};

#endif  // GRAPHICSINSTANCEITEM_H
