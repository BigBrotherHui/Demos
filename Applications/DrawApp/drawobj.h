#ifndef DRAWOBJ
#define DRAWOBJ

#include <qgraphicsitem.h>
#include <QCursor>
#include <QGraphicsItemGroup>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QList>
#include <QMimeData>
#include <QPainter>
#include <QStyleOptionGraphicsItem>
#include <QVector>
#include <QXmlStreamReader>
#include <QtCore/QObject>
#include <vector>
#include "graphicsitem.h"
#include "graphicspolygonitem.h"
#include "graphicsrectitem.h"
#include "layermanager.h"
#include "sizehandle.h"
class Group;
class ShapeMimeData : public QMimeData {
  Q_OBJECT
 public:
  ShapeMimeData(QList<QGraphicsItem *> items);
  ~ShapeMimeData();
  QList<QGraphicsItem *> items() const;

 private:
  QList<QGraphicsItem *> m_items;
};

// class GraphicsEllipseItem : public GraphicsRectItem {
// public:
//  GraphicsEllipseItem(const QRect &rect, QGraphicsItem *parent = 0);
//  QPainterPath shape() const;
//  void control(int dir, const QPointF &delta);
//  QRectF boundingRect() const;
//  QGraphicsItem *duplicate() const;
//  QString displayName() const { return tr("ellipse"); }
//  virtual bool loadFromXml(QXmlStreamReader *xml);
//  virtual bool saveToXml(QXmlStreamWriter *xml);

// protected:
//  void updatehandles();
//  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
//             QWidget *widget);
//  int m_startAngle;
//  int m_spanAngle;
//};

#endif  // DRAWOBJ
