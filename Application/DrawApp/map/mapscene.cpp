#include "mapscene.h"
#include <QPainter>
MapScene::MapScene(QObject *parent) : QGraphicsScene(parent) {}

void MapScene::drawBackground(QPainter *painter, const QRectF &rect) {
  QGraphicsScene::drawBackground(painter, rect);
  painter->fillRect(sceneRect(), Qt::white);
  painter->drawPixmap(rect, m_pix, rect);
}
