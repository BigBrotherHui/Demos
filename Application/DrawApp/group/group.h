#ifndef GROUP_H
#define GROUP_H
#include <QObject>
#include <QSet>
class GraphicsItem;
class Group : public QObject {
  Q_OBJECT
 public:
  Group(QObject *obj = nullptr);
  void addToGroup(GraphicsItem *item);
  void removeFromGroup(GraphicsItem *item);
  QSet<GraphicsItem *> items() { return m_items; }
  void selectAll(bool p);
  void addToGroup(QSet<GraphicsItem *> items);

 private:
  QSet<GraphicsItem *> m_items;
};

#endif  // GROUP_H
