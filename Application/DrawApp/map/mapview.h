#ifndef MAPVIEW_H
#define MAPVIEW_H

#include <QGraphicsView>
class BoxTool;
class MapScene;
class MapView : public QGraphicsView {
 public:
  MapView(QGraphicsView* map_view, QWidget* parent = nullptr);
 public slots:
  void updateImage();

 protected:
  void mouseMoveEvent(QMouseEvent* event) override;

 private:
  QGraphicsView* map_view;
  BoxTool* item;
  MapScene* m_scene;
};

#endif  // MAPVIEW_H
