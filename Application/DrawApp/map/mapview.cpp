#include "mapview.h"
#include <QDebug>
#include <QGraphicsItem>
#include <QPainter>
#include "boxtool.h"
#include "mapscene.h"
MapView::MapView(QGraphicsView *map_view, QWidget *parent)
    : QGraphicsView(parent), map_view(map_view) {
  setDragMode(QGraphicsView::NoDrag);
  m_scene = new MapScene;
  m_scene->setSceneRect(0, 0, 200, 200);
  setScene(m_scene);
  setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
  setRenderHint(QPainter::Antialiasing);
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
  targetR.setTop((sceneRect().height() - pos.y()) / sceneRect().height() *
                 map_view->sceneRect().height());
  //  targetR.setLeft(pos.x() / sceneRect().width() *
  //                  map_view->sceneRect().width());
  //  targetR.setTop((size().height() - pos.y() - item->rect().height()) /
  //                 sceneRect().height() * map_view->sceneRect().height());
  targetR.setWidth(800);
  targetR.setHeight(800);
  map_view->fitInView(targetR, Qt::KeepAspectRatio);
}
