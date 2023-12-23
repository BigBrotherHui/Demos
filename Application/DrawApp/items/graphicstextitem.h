#ifndef GRAPHICSTEXTITEM_H
#define GRAPHICSTEXTITEM_H

#include "graphicsrectitem.h"
class QGraphicsProxyWidget;
class GraphicsTextItem : public GraphicsRectItem {
 public:
  GraphicsTextItem(const QRect &rect, QGraphicsItem *parent = 0);
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
             QWidget *widget);
  virtual bool loadFromXml(QXmlStreamReader *xml);
  virtual bool saveToXml(QXmlStreamWriter *xml);
  QGraphicsProxyWidget *pMyProxy;
};

#endif  // GRAPHICSTEXTITEM_H
