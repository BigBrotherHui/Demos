#ifndef MAPSCENE_H
#define MAPSCENE_H

#include <QGraphicsScene>
class MapScene : public QGraphicsScene {
  Q_OBJECT
 public:
  MapScene(QObject* parent = nullptr);
  void drawBackground(QPainter* painter, const QRectF& rect) override;
 public slots:
  void updateImage(QPixmap pix) {
    m_pix = pix;
    update();
  }

 private:
  QPixmap m_pix;
};

#endif  // MAPSCENE_H
