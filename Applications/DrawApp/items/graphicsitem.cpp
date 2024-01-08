#include "graphicsitem.h"
#include <QDebug>
#include <QPainter>
#include <QStyleOptionGraphicsItem>
#include <QUuid>
GraphicsItem::GraphicsItem(QGraphicsItem *parent)
    : AbstractShapeType<QGraphicsItem>(parent) {
  /*
  QGraphicsDropShadowEffect *effect = new QGraphicsDropShadowEffect;
  effect->setBlurRadius(4);
  setGraphicsEffect(effect);
 */
  // handles
  setFlag(QGraphicsItem::ItemIsMovable, true);
  setFlag(QGraphicsItem::ItemIsSelectable, true);
  setFlag(QGraphicsItem::ItemSendsGeometryChanges, true);
  this->setAcceptHoverEvents(true);
  setObjectName(QUuid::createUuid().toString());
}

QPixmap GraphicsItem::image() {
  QPixmap pixmap(64, 64);
  pixmap.fill(Qt::transparent);
  QPainter painter(&pixmap);
  setPen(QPen(Qt::black));
  setBrush(Qt::white);
  QStyleOptionGraphicsItem *styleOption = new QStyleOptionGraphicsItem;
  //    painter.translate(m_localRect.center().x(),m_localRect.center().y());
  paint(&painter, styleOption);
  delete styleOption;
  return pixmap;
}

void GraphicsItem::addGroup(Group *g) {
  if (!m_groups.contains(g)) m_groups.insert(g);
}

void GraphicsItem::removeGroup(Group *g) {
  if (m_groups.contains(g)) m_groups.remove(g);
}

QSet<Group *> GraphicsItem::groups() { return m_groups; }

void GraphicsItem::mousePressEvent(QGraphicsSceneMouseEvent *event) {
  if (!iscreatehandles) {
    m_handles.reserve(Left);
    for (int i = LeftTop; i <= Left; ++i) {
      SizeHandleRect *shr = new SizeHandleRect(this, i);
      m_handles.push_back(shr);
    }
    iscreatehandles = 1;
  }
  return AbstractShapeType<QGraphicsItem>::mousePressEvent(event);
}

void GraphicsItem::mouseReleaseEvent(QGraphicsSceneMouseEvent *event) {
  return AbstractShapeType<QGraphicsItem>::mouseReleaseEvent(event);
}

void GraphicsItem::updatehandles() {
  const QRectF &geom = this->boundingRect();
  const Handles::iterator hend = m_handles.end();
  for (Handles::iterator it = m_handles.begin(); it != hend; ++it) {
    SizeHandleRect *hndl = *it;
    ;
    switch (hndl->dir()) {
      case LeftTop:
        hndl->move(geom.x(), geom.y());
        break;
      case Top:
        hndl->move(geom.x() + geom.width() / 2, geom.y());
        break;
      case RightTop:
        hndl->move(geom.x() + geom.width(), geom.y());
        break;
      case Right:
        hndl->move(geom.x() + geom.width(), geom.y() + geom.height() / 2);
        break;
      case RightBottom:
        hndl->move(geom.x() + geom.width(), geom.y() + geom.height());
        break;
      case Bottom:
        hndl->move(geom.x() + geom.width() / 2, geom.y() + geom.height());
        break;
      case LeftBottom:
        hndl->move(geom.x(), geom.y() + geom.height());
        break;
      case Left:
        hndl->move(geom.x(), geom.y() + geom.height() / 2);
        break;
      default:
        break;
    }
  }
}

void GraphicsItem::contextMenuEvent(QGraphicsSceneContextMenuEvent *event) {
  Q_UNUSED(event);
}

bool GraphicsItem::readBaseAttributes(QXmlStreamReader *xml) {
  qreal x = xml->attributes().value(tr("x")).toDouble();
  qreal y = xml->attributes().value(tr("y")).toDouble();
  m_width = xml->attributes().value("width").toDouble();
  m_height = xml->attributes().value("height").toDouble();
  setZValue(xml->attributes().value("z").toDouble());
  //  setRotation(xml->attributes().value("rotate").toDouble());
  setPos(x, y);
  return true;
}

bool GraphicsItem::writeBaseAttributes(QXmlStreamWriter *xml) {
  //  xml->writeAttribute(tr("rotate"), QString("%1").arg(rotation()));
  xml->writeAttribute(tr("x"), QString("%1").arg(pos().x()));
  xml->writeAttribute(tr("y"), QString("%1").arg(pos().y()));
  xml->writeAttribute(tr("z"), QString("%1").arg(zValue()));
  xml->writeAttribute(tr("width"), QString("%1").arg(m_width));
  xml->writeAttribute(tr("height"), QString("%1").arg(m_height));
  return true;
}

QVariant GraphicsItem::itemChange(QGraphicsItem::GraphicsItemChange change,
                                  const QVariant &value) {
  if (change == QGraphicsItem::ItemSelectedHasChanged) {
    QGraphicsItemGroup *g = dynamic_cast<QGraphicsItemGroup *>(parentItem());
    if (!g)
      setState(value.toBool() ? SelectionHandleActive : SelectionHandleOff);
    else {
      setSelected(false);
      return QVariant::fromValue<bool>(false);
    }
  }
  /*
  else if (change == ItemPositionChange && scene()) {
      // value is the new position.
      QPointF newPos = value.toPointF();
      QRectF rect = scene()->sceneRect();
      if (!rect.contains(newPos)) {
          // Keep the item inside the scene rect.
          newPos.setX(qMin(rect.right()-boundingRect().width()/2,
  qMax(newPos.x(), rect.left()+boundingRect().width()/2)));
          newPos.setY(qMin(rect.bottom()-boundingRect().height()/2,
  qMax(newPos.y(), rect.top()+boundingRect().height()/2))); return newPos;
      }
  }
  */
  return QGraphicsItem::itemChange(change, value);
}
