#include "graphicsrectitem.h"
#include <QDebug>
#include <QPainter>
#include <QStyleOptionGraphicsItem>
#include "drawscene.h"
#include "util.h"
#include "globalsignalinstance.h"
GraphicsRectItem::GraphicsRectItem(const QRect &rect, bool isRound,
                                   QGraphicsItem *parent)
    : GraphicsItem(parent),
      m_isRound(isRound),
      m_fRatioX(1 / 10.0),
      m_fRatioY(1 / 3.0) {
  m_itemtype = rectangle;
  m_width = rect.width();
  m_height = rect.height();
  m_initialRect = rect;
  m_localRect = m_initialRect;
  m_localRect = rect;
  m_originPoint = QPointF(0, 0);
  updatehandles();
  connect(this,&GraphicsRectItem::refreshItem,GlobalSignalInstance::instance(),&GlobalSignalInstance::signal_wrapper_refreshItem);
}

QRectF GraphicsRectItem::boundingRect() const { return m_localRect; }

QPainterPath GraphicsRectItem::shape() const {
  QPainterPath path;
  double rx, ry;
  if (m_fRatioX <= 0)
    rx = 0;
  else {
    rx = m_width * m_fRatioX + 0.5;
  }
  if (m_fRatioY <= 0)
    ry = 0;
  else
    ry = m_height * m_fRatioY + 0.5;
  if (m_isRound)
    path.addRoundedRect(rect(), rx, ry);
  else
    path.addRect(rect());
  return path;
}

void GraphicsRectItem::control(int dir, const QPointF &delta) {
  QPointF local = mapFromParent(delta);
  switch (dir) {
    case 9: {
      QRectF delta1 = rect();
      int y = local.y();
      if (y > delta1.center().y()) y = delta1.center().y();
      if (y < delta1.top()) y = delta1.top();
      int H = delta1.height();
      if (H == 0) H = 1;
      m_fRatioY = std::abs(((float)(delta1.top() - y))) / H;
    } break;
    case 10: {
      QRectF delta1 = rect();
      int x = local.x();
      if (x < delta1.center().x()) x = delta1.center().x();
      if (x > delta1.right()) x = delta1.right();
      int W = delta1.width();
      if (W == 0) W = 1;
      m_fRatioX = std::abs(((float)(delta1.right() - x))) / W;
      break;
    }
    case 11: {
      //        setTransform(transform().translate(-local.x(),-local.y()));
      //        setTransformOriginPoint(local.x(),local.y());
      //        setTransform(transform().translate(local.x(),local.y()));
      m_originPoint = local;
    } break;
    default:
      break;
  }
  prepareGeometryChange();
  updatehandles();
}

void GraphicsRectItem::stretch(int handle, double sx, double sy,
                               const QPointF &origin) {
  QTransform trans;
  switch (handle) {
    case Right:
    case Left:
      sy = 1;
      break;
    case Top:
    case Bottom:
      sx = 1;
      break;
    default:
      break;
  }

  opposite_ = origin;

  trans.translate(origin.x(), origin.y());
  trans.scale(sx, sy);
  trans.translate(-origin.x(), -origin.y());

  prepareGeometryChange();
  m_localRect = trans.mapRect(m_initialRect);
  m_width = m_localRect.width();
  m_height = m_localRect.height();
  updatehandles();
}

void GraphicsRectItem::updateCoordinate() {
  QPointF pt1, pt2, delta;
  pt1 = mapToScene(transformOriginPoint());
  pt2 = mapToScene(m_localRect.center());
  pt1 = pt1.toPoint() / DrawScene::m_mingridsize.width() *
        DrawScene::m_mingridsize.width();
  pt2 = pt2.toPoint() / DrawScene::m_mingridsize.height() *
        DrawScene::m_mingridsize.height();
  delta = pt1 - pt2;
  delta = delta.toPoint() / DrawScene::m_mingridsize.height() *
          DrawScene::m_mingridsize.height();
  if (!parentItem()) {
//    prepareGeometryChange();
//    int x = -m_width / 2;
//    int y = -m_height / 2;
//    x = x / DrawScene::m_mingridsize.width() * DrawScene::m_mingridsize.width();
//    y = y / DrawScene::m_mingridsize.height() *DrawScene::m_mingridsize.height();
    m_width = qRound(m_width) / DrawScene::m_mingridsize.width() *
              DrawScene::m_mingridsize.width();
    m_height = qRound(m_height) / DrawScene::m_mingridsize.height() *
               DrawScene::m_mingridsize.height();
//    m_localRect = QRect(x, y, m_width, m_height);
//    // TODO:这个地方如果新添加polygon或polyline会导致第一个矩形反向绘制
//    // TODO:刻度尺y轴值需要取反
//    setTransform(transform().translate(delta.x(), delta.y()));
//    setTransformOriginPoint(m_localRect.center());
//    moveBy(-delta.x(), -delta.y());
//    setTransform(transform().translate(-delta.x(), -delta.y()));
//    opposite_ = QPointF(0, 0);
    updatehandles();
    emit refreshItem(this);
  }
  m_initialRect = m_localRect;
}

void GraphicsRectItem::move(const QPointF &point) {
  moveBy(point.x(), point.y());
}

QGraphicsItem *GraphicsRectItem::duplicate() const {
  GraphicsRectItem *item = new GraphicsRectItem(rect().toRect(), m_isRound);
  item->createHandles();
  item->m_width = width();
  item->m_height = height();
  item->setPos(pos().x(), pos().y());
  item->setPen(pen());
  item->setBrush(brush());
  item->setTransform(transform());
  item->setTransformOriginPoint(transformOriginPoint());
  item->setRotation(rotation());
  item->setScale(scale());
  item->setZValue(zValue());
  item->m_fRatioY = m_fRatioY;
  item->m_fRatioX = m_fRatioX;
  item->updateCoordinate();
  return item;
}

bool GraphicsRectItem::loadFromXml(QXmlStreamReader *xml) {
  m_isRound = (xml->name() == tr("roundrect"));
  if (m_isRound) {
    m_fRatioX = xml->attributes().value(tr("rx")).toDouble();
    m_fRatioY = xml->attributes().value(tr("ry")).toDouble();
  }
  readBaseAttributes(xml);
  updateCoordinate();
  xml->skipCurrentElement();
  return true;
}

bool GraphicsRectItem::saveToXml(QXmlStreamWriter *xml) {
  if (m_isRound) {
    xml->writeStartElement(tr("roundrect"));
    xml->writeAttribute(tr("rx"), QString("%1").arg(m_fRatioX));
    xml->writeAttribute(tr("ry"), QString("%1").arg(m_fRatioY));
  } else
    xml->writeStartElement(tr("rect"));

  writeBaseAttributes(xml);
  xml->writeEndElement();
  return true;
}

void GraphicsRectItem::updatehandles() {
  const QRectF &geom = this->boundingRect();
  GraphicsItem::updatehandles();
  if (m_isRound) {
    m_handles[8]->move(geom.right(), geom.top() + geom.height() * m_fRatioY);
    m_handles[9]->move(geom.right() - geom.width() * m_fRatioX, geom.top());
    // m_handles[10]->move(m_originPoint.x(),m_originPoint.y());
  }
}
void GraphicsRectItem::paint(QPainter *painter,
                             const QStyleOptionGraphicsItem *option,
                             QWidget *widget) {
  painter->setPen(pen());
  painter->setBrush(brush());
  double rx, ry;
  if (m_fRatioX <= 0)
    rx = 0;
  else {
    rx = m_width * m_fRatioX + 0.5;
  }
  if (m_fRatioY <= 0)
    ry = 0;
  else
    ry = m_height * m_fRatioY + 0.5;
  if (m_isRound)
    painter->drawRoundedRect(rect(), rx, ry);
  else
    painter->drawRect(rect().toRect());
  painter->setPen(Qt::blue);
  if (option->state & QStyle::State_Selected)
    qt_graphicsItem_highlightSelected(this, painter, option);
}
