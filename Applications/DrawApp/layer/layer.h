#ifndef LAYER_H
#define LAYER_H

#include <QBrush>
#include <QObject>
#include <QPen>
#include <vector>
class QXmlStreamReader;
class QXmlStreamWriter;
class DrawScene;
class GraphicsItem;
class QGraphicsItem;
class Layer : public QObject {
  Q_OBJECT
  Q_PROPERTY(QColor m_pen READ penColor WRITE setPenColor)
  Q_PROPERTY(QColor m_brush READ brushColor WRITE setBrushColor)
  Q_PROPERTY(double m_z READ GetZValue WRITE SetZValue)
 public:
  ~Layer();
  void SetZValue(int z);
  void AddItem(GraphicsItem *item);
  void RemoveItem(GraphicsItem *item);
  int GetZValue() { return m_z; }
  QList<GraphicsItem *> GetAllItems();
  void SetVisible(bool b);
  void SetMovable(bool e);
  void SetSelectable(bool e);
  void SetScene(DrawScene *s) { m_scene = s; }
  DrawScene *Scene() { return m_scene; }
  void selectAll();
  QColor brushColor() { return m_brush.color(); }
  void setBrushColor(QColor c);
  QColor penColor() { return m_pen.color(); }
  void setPenColor(QColor c);
  Layer *loadFromXml(QXmlStreamReader *xml);
  bool saveToXml(QXmlStreamWriter *xml);
  GraphicsItem *selectItemByObjectName(QString name);

 protected:
  explicit Layer(int z, DrawScene *s);

 private:
  int m_z;
  QList<GraphicsItem *> m_items;
  QBrush m_brush{Qt::red};
  QPen m_pen{Qt::green};
  DrawScene *m_scene{nullptr};
  friend class LayerManager;
};

#endif  // LAYER_H
