#ifndef UTIL_H
#define UTIL_H
#include <QGraphicsItem>
#include <QPainter>
#include <QPainterPath>
#include <QStyleOptionGraphicsItem>
QPainterPath qt_graphicsItem_shapeFromPath(const QPainterPath &path,
                                           const QPen &pen);

void qt_graphicsItem_highlightSelected(QGraphicsItem *item, QPainter *painter,
                                       const QStyleOptionGraphicsItem *option);

#endif  // UTIL_H
