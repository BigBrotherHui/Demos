#ifndef LAYERMANAGER_H
#define LAYERMANAGER_H

#include <QObject>
class Layer;
class DrawScene;
class GraphicsItem;
class QTreeWidget;
class QTreeWidgetItem;
class QGraphicsItem;
class QXmlStreamReader;
class QListWidget;
class QTreeWidget;

class LayerManager : public QObject {
  Q_OBJECT
 public:
//  static LayerManager *GetInstance();
    LayerManager(DrawScene *scene,QObject *parent = nullptr);
    ~LayerManager();
  void SetItemWidget(QListWidget *w);
  Layer *AddLayer();
  void AddLayer(Layer *l);
  void RemoveLayer(Layer *);
  void SetCurrentLayer(Layer *l) {
    if (m_curlayer != l) m_curlayer = l;
  }
  bool isCurrentLayerValid();
  Layer *GetCurrentLayer() { return m_curlayer; }
  QList<Layer *> GetAllLayers() { return m_layers; }
  void AddItem(GraphicsItem *i);
  void AddToSpecifiedLayer(int z, GraphicsItem *item);
  static Layer *loadFromXml(QXmlStreamReader *r, DrawScene *scene);
  void refreshItemWidget();
  void selectItemByType(QString t);
  QTreeWidget *layerwidget{nullptr};
  void refreshLayerWidget();
 signals:
  void layerChanged();

 protected:
  Layer *GetLayer(int rowNo);
  void slot_cellClicked(QTreeWidgetItem *item, int);
  void SetLayerWidget(QTreeWidget *w);

 private:
  QList<Layer *> m_layers;
  Layer *m_curlayer{nullptr};
  int maxlayer = -1;
  DrawScene *m_scene{nullptr};
  QListWidget *itemwidget{nullptr};
  std::map<QTreeWidgetItem *, Layer *> layermap;
};

#endif  // LAYERMANAGER_H
