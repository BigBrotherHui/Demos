#ifndef BOXTOOL_H
#define BOXTOOL_H

#include <QGraphicsRectItem>
class BoxTool : public QGraphicsRectItem {
 public:
  BoxTool();

 protected:
  void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
             QWidget* widget) override;
  QVariant itemChange(GraphicsItemChange change,
                      const QVariant& value) override;
};

#endif  // BOXTOOL_H
