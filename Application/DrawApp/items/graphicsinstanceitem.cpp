#include "graphicsinstanceitem.h"
#include <QPainter>
GraphicsInstanceItem::GraphicsInstanceItem(const QRect &rect,
                                           QGraphicsItem *parent)
    : GraphicsRectItem(rect, parent) {}

QRectF GraphicsInstanceItem::boundingRect() const { return rect(); }

QPainterPath GraphicsInstanceItem::shape() const {
  QPainterPath pa;
  pa.addRect(rect());
  return pa;
}

void GraphicsInstanceItem::paint(QPainter *painter,
                                 const QStyleOptionGraphicsItem *option,
                                 QWidget *widget) {
  pixmap = pixmap.scaled(rect().size().toSize(), Qt::KeepAspectRatio);
  painter->drawPixmap(rect(), pixmap, rect());
}
