#include "util.h"
QPainterPath qt_graphicsItem_shapeFromPath(const QPainterPath &path,
                                           const QPen &pen) {
  // We unfortunately need this hack as QPainterPathStroker will set a width
  // of 1.0 if we pass a value of 0.0 to QPainterPathStroker::setWidth()
  const qreal penWidthZero = qreal(0.00000001);

  if (path == QPainterPath() || pen == Qt::NoPen) return path;
  QPainterPathStroker ps;
  ps.setCapStyle(pen.capStyle());
  if (pen.widthF() <= 0.0)
    ps.setWidth(penWidthZero);
  else
    ps.setWidth(pen.widthF());
  ps.setJoinStyle(pen.joinStyle());
  ps.setMiterLimit(pen.miterLimit());
  QPainterPath p = ps.createStroke(path);
  p.addPath(path);
  return p;
}

void qt_graphicsItem_highlightSelected(QGraphicsItem *item, QPainter *painter,
                                       const QStyleOptionGraphicsItem *option) {
  const QRectF murect = painter->transform().mapRect(QRectF(0, 0, 1, 1));
  if (qFuzzyIsNull(qMax(murect.width(), murect.height()))) return;

  const QRectF mbrect = painter->transform().mapRect(item->boundingRect());
  if (qMin(mbrect.width(), mbrect.height()) < qreal(1.0)) return;

  qreal itemPenWidth;
  switch (item->type()) {
      //    case QGraphicsEllipseItem::Type:
      //      itemPenWidth = static_cast<QGraphicsEllipseItem
      //      *>(item)->pen().widthF(); break;
    case QGraphicsPathItem::Type:
      itemPenWidth = static_cast<QGraphicsPathItem *>(item)->pen().widthF();
      break;
    case QGraphicsPolygonItem::Type:
      itemPenWidth = static_cast<QGraphicsPolygonItem *>(item)->pen().widthF();
      break;
    case QGraphicsRectItem::Type:
      itemPenWidth = static_cast<QGraphicsRectItem *>(item)->pen().widthF();
      break;
    case QGraphicsSimpleTextItem::Type:
      itemPenWidth =
          static_cast<QGraphicsSimpleTextItem *>(item)->pen().widthF();
      break;
    case QGraphicsLineItem::Type:
      itemPenWidth = static_cast<QGraphicsLineItem *>(item)->pen().widthF();
      break;
    default:
      itemPenWidth = 1.0;
  }
  const qreal pad = itemPenWidth / 2;

  const qreal penWidth = 0;  // cosmetic pen

  const QColor fgcolor = option->palette.windowText().color();
  const QColor bgcolor(  // ensure good contrast against fgcolor
      fgcolor.red() > 127 ? 0 : 255, fgcolor.green() > 127 ? 0 : 255,
      fgcolor.blue() > 127 ? 0 : 255);

  painter->setPen(QPen(bgcolor, penWidth, Qt::SolidLine));
  painter->setBrush(Qt::NoBrush);
  painter->drawRect(item->boundingRect().adjusted(-pad, -pad, pad, pad));

  painter->setPen(QPen(QColor("lightskyblue"), 0, Qt::SolidLine));
  painter->setBrush(Qt::NoBrush);
  painter->drawRect(item->boundingRect().adjusted(-pad, -pad, pad, pad));
}
