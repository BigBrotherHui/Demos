#include "drawobj.h"
#include <QCursor>
#include <QDateTime>
#include <QDebug>
#include <QGraphicsEffect>
#include <QGraphicsProxyWidget>
#include <QGraphicsTransform>
#include <QLinearGradient>
#include <QMatrix4x4>
#include <QPainter>
#include <QPlainTextEdit>
#include <QRadialGradient>
#include <QStyle>
#include <QStyleOptionGraphicsItem>
#include <QUuid>
#include <cmath>
#include "drawscene.h"

#include "group.h"
#define M_PI 3.1415926
ShapeMimeData::ShapeMimeData(QList<QGraphicsItem *> items) {
  foreach (QGraphicsItem *item, items) {
    AbstractShape *sp = qgraphicsitem_cast<AbstractShape *>(item);
    m_items.append(sp->duplicate());
  }
}
ShapeMimeData::~ShapeMimeData() {
  foreach (QGraphicsItem *item, m_items) { delete item; }
  m_items.clear();
}

QList<QGraphicsItem *> ShapeMimeData::items() const { return m_items; }

// void qt_graphicsItem_highlightSelected(QGraphicsItem *item, QPainter
// *painter,
//                                       const QStyleOptionGraphicsItem *option)
//                                       {
//  const QRectF murect = painter->transform().mapRect(QRectF(0, 0, 1, 1));
//  if (qFuzzyIsNull(qMax(murect.width(), murect.height()))) return;

//  const QRectF mbrect = painter->transform().mapRect(item->boundingRect());
//  if (qMin(mbrect.width(), mbrect.height()) < qreal(1.0)) return;

//  qreal itemPenWidth;
//  switch (item->type()) {
//      //    case QGraphicsEllipseItem::Type:
//      //      itemPenWidth = static_cast<QGraphicsEllipseItem
//      //      *>(item)->pen().widthF(); break;
//    case QGraphicsPathItem::Type:
//      itemPenWidth = static_cast<QGraphicsPathItem *>(item)->pen().widthF();
//      break;
//    case QGraphicsPolygonItem::Type:
//      itemPenWidth = static_cast<QGraphicsPolygonItem
//      *>(item)->pen().widthF(); break;
//    case QGraphicsRectItem::Type:
//      itemPenWidth = static_cast<QGraphicsRectItem *>(item)->pen().widthF();
//      break;
//    case QGraphicsSimpleTextItem::Type:
//      itemPenWidth =
//          static_cast<QGraphicsSimpleTextItem *>(item)->pen().widthF();
//      break;
//    case QGraphicsLineItem::Type:
//      itemPenWidth = static_cast<QGraphicsLineItem *>(item)->pen().widthF();
//      break;
//    default:
//      itemPenWidth = 1.0;
//  }
//  const qreal pad = itemPenWidth / 2;

//  const qreal penWidth = 0;  // cosmetic pen

//  const QColor fgcolor = option->palette.windowText().color();
//  const QColor bgcolor(  // ensure good contrast against fgcolor
//      fgcolor.red() > 127 ? 0 : 255, fgcolor.green() > 127 ? 0 : 255,
//      fgcolor.blue() > 127 ? 0 : 255);

//  painter->setPen(QPen(bgcolor, penWidth, Qt::SolidLine));
//  painter->setBrush(Qt::NoBrush);
//  painter->drawRect(item->boundingRect().adjusted(-pad, -pad, pad, pad));

//  painter->setPen(QPen(QColor("lightskyblue"), 0, Qt::SolidLine));
//  painter->setBrush(Qt::NoBrush);
//  painter->drawRect(item->boundingRect().adjusted(-pad, -pad, pad, pad));
//}

static QRectF RecalcBounds(const QPolygonF &pts) {
  QRectF bounds(pts[0], QSize(0, 0));
  for (int i = 1; i < pts.count(); ++i) {
    if (pts[i].x() < bounds.left()) bounds.setLeft(pts[i].x());
    if (pts[i].x() > bounds.right()) bounds.setRight(pts[i].x());
    if (pts[i].y() < bounds.top()) bounds.setTop(pts[i].y());
    if (pts[i].y() > bounds.bottom()) bounds.setBottom(pts[i].y());
  }
  bounds = bounds.normalized();
  return bounds;
}

// GraphicsEllipseItem::GraphicsEllipseItem(const QRect &rect,
//                                         QGraphicsItem *parent)
//    : GraphicsRectItem(rect, parent) {
//  m_startAngle = 40;
//  m_spanAngle = 400;
//  SizeHandleRect *shr = new SizeHandleRect(this, 9, true);
//  m_handles.push_back(shr);
//  shr = new SizeHandleRect(this, 10, true);
//  m_handles.push_back(shr);
//  updatehandles();
//}

// QPainterPath GraphicsEllipseItem::shape() const {
//  QPainterPath path;
//  int startAngle = m_startAngle <= m_spanAngle ? m_startAngle : m_spanAngle;
//  int endAngle = m_startAngle >= m_spanAngle ? m_startAngle : m_spanAngle;
//  if (endAngle - startAngle > 360) endAngle = startAngle + 360;

//  if (m_localRect.isNull()) return path;
//  if ((endAngle - startAngle) % 360 != 0) {
//    path.moveTo(m_localRect.center());
//    path.arcTo(m_localRect, startAngle, endAngle - startAngle);
//  } else {
//    path.addEllipse(m_localRect);
//  }
//  path.closeSubpath();
//  return path;
//}

// void GraphicsEllipseItem::control(int dir, const QPointF &delta) {
//  QPointF local = mapFromScene(delta);

//  switch (dir) {
//    case 9: {
//      qreal len_y = local.y() - m_localRect.center().y();
//      qreal len_x = local.x() - m_localRect.center().x();
//      m_startAngle = -atan2(len_y, len_x) * 180 / M_PI;
//    } break;
//    case 10: {
//      qreal len_y = local.y() - m_localRect.center().y();
//      qreal len_x = local.x() - m_localRect.center().x();
//      m_spanAngle = -atan2(len_y, len_x) * 180 / M_PI;
//      break;
//    }
//    default:
//      break;
//  }
//  prepareGeometryChange();
//  if (m_startAngle > m_spanAngle) m_startAngle -= 360;
//  if (m_spanAngle < m_startAngle) {
//    qreal tmp = m_spanAngle;
//    m_spanAngle = m_startAngle;
//    m_startAngle = tmp;
//  }

//  if (qAbs(m_spanAngle - m_startAngle) > 360) {
//    m_startAngle = 40;
//    m_spanAngle = 400;
//  }
//  updatehandles();
//}

// QRectF GraphicsEllipseItem::boundingRect() const {
//  return shape().controlPointRect();
//}

// QGraphicsItem *GraphicsEllipseItem::duplicate() const {
//  GraphicsEllipseItem *item = new GraphicsEllipseItem(m_localRect.toRect());
//  item->m_width = width();
//  item->m_height = height();
//  item->m_startAngle = m_startAngle;
//  item->m_spanAngle = m_spanAngle;

//  item->setPos(pos().x(), pos().y());
//  item->setPen(pen());
//  item->setBrush(brush());
//  item->setTransform(transform());
//  item->setTransformOriginPoint(transformOriginPoint());
//  item->setRotation(rotation());
//  item->setScale(scale());
//  item->setZValue(zValue() + 0.1);
//  item->updateCoordinate();
//  return item;
//}

// bool GraphicsEllipseItem::loadFromXml(QXmlStreamReader *xml) {
//  m_startAngle = xml->attributes().value("startAngle").toInt();
//  m_spanAngle = xml->attributes().value("spanAngle").toInt();
//  readBaseAttributes(xml);
//  xml->skipCurrentElement();
//  updateCoordinate();
//  return true;
//}

// bool GraphicsEllipseItem::saveToXml(QXmlStreamWriter *xml) {
//  xml->writeStartElement(tr("ellipse"));
//  xml->writeAttribute("startAngle", QString("%1").arg(m_startAngle));
//  xml->writeAttribute("spanAngle", QString("%1").arg(m_spanAngle));

//  writeBaseAttributes(xml);
//  xml->writeEndElement();
//  return true;
//}

// void GraphicsEllipseItem::updatehandles() {
//  GraphicsItem::updatehandles();
//  QRectF local = QRectF(-m_width / 2, -m_height / 2, m_width, m_height);
//  QPointF delta = local.center() - m_localRect.center();

//  qreal x = (m_width / 2) * cos(-m_startAngle * M_PI / 180);
//  qreal y = (m_height / 2) * sin(-m_startAngle * M_PI / 180);

//  m_handles.at(8)->move(x - delta.x(), y - delta.y());
//  x = (m_width / 2) * cos(-m_spanAngle * M_PI / 180);
//  y = (m_height / 2) * sin(-m_spanAngle * M_PI / 180);
//  m_handles.at(9)->move(x - delta.x(), y - delta.y());
//}

// void GraphicsEllipseItem::paint(QPainter *painter,
//                                const QStyleOptionGraphicsItem *option,
//                                QWidget *widget) {
//  Q_UNUSED(option);
//  Q_UNUSED(widget);
//  QColor c = brushColor();
//  QRectF rc = m_localRect;

//  qreal radius = qMax(rc.width(), rc.height());

//  /*
//      QConicalGradient  result(rc.center(),-45);

//      QColor niceBlue(150, 150, 200);
//      result.setColorAt(0.0, c.dark(200));
//      result.setColorAt(0.2, niceBlue);
//      result.setColorAt(0.5, c.light(120));
//      result.setColorAt(1.0, c.dark(200));
//  */
//  painter->setPen(pen());
//  QBrush b(c);
//  painter->setBrush(b);

//  int startAngle = m_startAngle <= m_spanAngle ? m_startAngle : m_spanAngle;
//  int endAngle = m_startAngle >= m_spanAngle ? m_startAngle : m_spanAngle;
//  if (endAngle - startAngle > 360) endAngle = startAngle + 360;

//  if (qAbs(endAngle - startAngle) % (360) == 0)
//    painter->drawEllipse(m_localRect);
//  else
//    painter->drawPie(m_localRect, startAngle * 16,
//                     (endAngle - startAngle) * 16);

//  if (option->state & QStyle::State_Selected)
//    qt_graphicsItem_highlightSelected(this, painter, option);
//}
