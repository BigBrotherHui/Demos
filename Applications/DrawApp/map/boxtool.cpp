#include "boxtool.h"
#include <QGraphicsScene>
#include <QPainter>
#include <QPen>
#include <QStyleOptionGraphicsItem>
BoxTool::BoxTool() { setFlag(QGraphicsItem::ItemIsMovable, true); }

void BoxTool::paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
                    QWidget *widget) {
  QStyleOptionGraphicsItem op;
  op.initFrom(widget);
  if (option->state & QStyle::State_Selected) op.state = QStyle::State_None;

  QGraphicsRectItem::paint(painter, &op, widget);

  if (option->state & QStyle::State_Selected) {
    qreal itemPenWidth = pen().widthF();
    const qreal pad = itemPenWidth;
    const qreal penWidth = .2;
    QColor color = QColor(Qt::green);
    painter->setPen(QPen(color, penWidth, Qt::SolidLine));
    painter->setBrush(Qt::NoBrush);
    painter->drawRect(boundingRect());
  }
}

QVariant BoxTool::itemChange(QGraphicsItem::GraphicsItemChange change,
                             const QVariant &value) {
  if (change == ItemPositionChange && scene()) {
    QPointF newPos = value.toPointF();
    QRectF sceneRect = scene()->sceneRect();  //水平移动
    sceneRect.adjust(0, 0, -rect().width(), -rect().height());
    if (!sceneRect.contains(newPos)) {
      newPos.setX(qMin(sceneRect.right(), qMax(newPos.x(), sceneRect.left())));
      newPos.setY(qMin(sceneRect.bottom(), qMax(newPos.y(), sceneRect.top())));
      return newPos;
    }
  }
  return QGraphicsRectItem::itemChange(change, value);
}
