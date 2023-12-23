#include "mapview.h"
#include <QDebug>
#include <QGraphicsItem>
#include <QPainter>
#include "boxtool.h"
#include "mapscene.h"
#include <QTransform>
MapView::MapView(QGraphicsView *map_view, QWidget *parent)
    : QGraphicsView(parent), map_view(map_view) {
  setDragMode(QGraphicsView::NoDrag);
  m_scene = new MapScene;
  m_scene->setSceneRect(0, 0, 200, 200);
  setScene(m_scene);
  setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
  setOptimizationFlags(QGraphicsView::DontSavePainterState);
  setViewportUpdateMode(QGraphicsView::SmartViewportUpdate);
  setRenderHint(QPainter::Antialiasing);
  setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  setFixedSize(202, 202);
  setAlignment(Qt::AlignLeft | Qt::AlignTop);
  setSceneRect(m_scene->sceneRect());

  item = new BoxTool;
  item->setPen(QPen(QColor(0, 0, 0), 2));
  item->setRect(0, 0, 40, 40);
  item->setPos(80, 80);
  item->setFlags(QGraphicsItem::ItemIsSelectable |
                 QGraphicsItem::ItemIsMovable |
                 QGraphicsItem::ItemSendsGeometryChanges);
  m_scene->addItem(item);
}

void MapView::updateImage() {
  QPixmap pix(size());
  QPainter painter(&pix);
  painter.setRenderHint(QPainter::Antialiasing);
//  item->scale(map_view->matrix().m11(),map_view->matrix().m11());
  map_view->render(&painter, QRect(0, 0, 400, 400), QRect(0, 0, 3200, 3200));
  m_scene->updateImage(pix);
}

void MapView::mouseMoveEvent(QMouseEvent *event) {
  QGraphicsView::mouseMoveEvent(event);
  if (!item || !item->isUnderMouse()) return;
  QPointF pos = mapToScene(item->pos().toPoint());
  QRect targetR;
  targetR.setLeft(pos.x() / sceneRect().width() *
                  map_view->sceneRect().width());
  targetR.setTop(pos.y() / sceneRect().height() *
                 map_view->sceneRect().height());
  //  targetR.setLeft(pos.x() / sceneRect().width() *
  //                  map_view->sceneRect().width());
  //  targetR.setTop((size().height() - pos.y() - item->rect().height()) /
  //                 sceneRect().height() * map_view->sceneRect().height());
  targetR.setWidth(5500);
  targetR.setHeight(2500);
  map_view->fitInView(targetR, Qt::KeepAspectRatio);
}
