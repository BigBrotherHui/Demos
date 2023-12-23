#ifndef DRAWSCENE
#define DRAWSCENE

#include <QGraphicsScene>
#include "drawobj.h"
#include "drawtool.h"

QT_BEGIN_NAMESPACE
class QGraphicsSceneMouseEvent;
class QMenu;
class QPointF;
class QGraphicsLineItem;
class QFont;
class QGraphicsTextItem;
class QColor;
class QKeyEvent;
class QPainter;
QT_END_NAMESPACE

enum AlignType {
  UP_ALIGN = 0,
  HORZ_ALIGN,
  VERT_ALIGN,
  DOWN_ALIGN,
  LEFT_ALIGN,
  RIGHT_ALIGN,
  CENTER_ALIGN,
  HORZEVEN_ALIGN,
  VERTEVEN_ALIGN,
  WIDTH_ALIGN,
  HEIGHT_ALIGN,
  ALL_ALIGN
};

class GridTool {
 public:
  GridTool(const QSize &grid = QSize(800, 800));
  void paintGrid(QPainter *painter, const QRect &rect);

 protected:
  QSize m_sizeGrid;
};

class GraphicsItemGroup;

class DrawScene : public QGraphicsScene {
  Q_OBJECT
 public:
  explicit DrawScene(QObject *parent = 0);
  ~DrawScene();
  void setView(QGraphicsView *view) { m_view = view; }
  QGraphicsView *view() { return m_view; }
  void align(AlignType alignType);
  void mouseEvent(QGraphicsSceneMouseEvent *mouseEvent);
  GraphicsItemGroup *createGroup(const QList<QGraphicsItem *> &items,
                                 bool isAdd = true);
  void destroyGroup(QGraphicsItemGroup *group);
  void SetMinGridSize(int s);
  void SetMinGridSize(QSize s);
  QSize gridSizeSpace() { return m_mingridsize; }
  static QSize m_mingridsize;
  bool getAllCollidingItems(QSet<GraphicsItem *> &all, GraphicsItem *item);
 signals:
  void itemMoved(QGraphicsItem *item, const QPointF &oldPosition);
  void itemRotate(QGraphicsItem *item, const qreal oldAngle);
  void itemAdded(QGraphicsItem *item);
  void itemResize(QGraphicsItem *item, int handle, const QPointF &scale);
  void itemControl(QGraphicsItem *item, int handle, const QPointF &newPos,
                   const QPointF &lastPos_);
  void itemRemoved();

 protected:
  void drawBackground(QPainter *painter, const QRectF &rect) Q_DECL_OVERRIDE;
  void mousePressEvent(QGraphicsSceneMouseEvent *mouseEvent) Q_DECL_OVERRIDE;
  void mouseMoveEvent(QGraphicsSceneMouseEvent *mouseEvent) Q_DECL_OVERRIDE;
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *mouseEvent) Q_DECL_OVERRIDE;
  void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *mouseEvet)
      Q_DECL_OVERRIDE;
  void focusInEvent(QFocusEvent *event) override;
  void focusOutEvent(QFocusEvent *event) override;
  void keyPressEvent(QKeyEvent *e) Q_DECL_OVERRIDE;
  void keyReleaseEvent(QKeyEvent *e) Q_DECL_OVERRIDE;
  QGraphicsView *m_view;
  qreal m_dx;
  qreal m_dy;
  bool m_moved;
  GridTool *m_grid;
  friend class Layer;
};

#endif  // DRAWSCENE
