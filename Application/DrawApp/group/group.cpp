#include "group.h"
#include <QUuid>
#include "drawobj.h"
Group::Group(QObject *obj) : QObject(obj) {
  setObjectName(QUuid::createUuid().toString());
}

void Group::addToGroup(GraphicsItem *item) {
  if (!m_items.contains(item)) {
    m_items.insert(item);
    item->addGroup(this);
  }
}

void Group::removeFromGroup(GraphicsItem *item) {
  if (m_items.contains(item)) {
    m_items.remove(item);
    item->removeGroup(this);
  }
  if (m_items.size() == 0) {
    delete this;
  }
}

void Group::selectAll(bool p) {
  for (auto i : m_items) {
    i->setSelected(p);
  }
}

void Group::addToGroup(QSet<GraphicsItem *> items) {
  m_items.unite(items);
  for (auto i : items) i->addGroup(this);
}
