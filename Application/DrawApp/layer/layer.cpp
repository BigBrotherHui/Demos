#include "layer.h"
#include <QDebug>
#include <QUuid>
#include "drawobj.h"
#include "drawscene.h"
Layer::Layer(int z, DrawScene *s) : m_z(z), m_scene(s) {
  setObjectName(QUuid::createUuid().toString());
}

Layer::~Layer() {
  for (int i = 0; i < m_items.size(); i++) {
    delete m_items[i];
  }
  m_items.clear();
}
void Layer::SetZValue(int z) {
  m_z = z;
  for (int i = 0; i < m_items.size(); i++) {
    m_items[i]->setZValue(m_z);
  }
}

void Layer::AddItem(GraphicsItem *item) {
  if (!item || m_items.contains(item)) return;
  m_items.push_back(item);
  item->setZValue(m_z);
  item->setBrush(m_brush);
  item->setPen(m_pen);
  m_scene->addItem(item);
  qDebug() << "Layer " << m_z << ":item size" << m_items.size();
}

void Layer::RemoveItem(GraphicsItem *item) {
  if (!item || !m_items.contains(item)) return;
  m_items.removeOne(item);
  if (m_scene) m_scene->removeItem(item);
  delete item;
  item = nullptr;
}

QList<GraphicsItem *> Layer::GetAllItems() {
  for (int i = 0; i < m_items.size(); i++) {
    if (!m_scene->items().contains(m_items[i])) {
      m_items.removeOne(m_items[i]);
    }
  }
  return m_items;
}

void Layer::SetVisible(bool b) {
  //  for (int i = 0; i < m_items.size(); i++) m_items[i]->SetVisible(b);
}

void Layer::SetMovable(bool e) {
  //  for (int i = 0; i < m_items.size(); i++) m_items[i]->SetMovable(e);
}

void Layer::SetSelectable(bool e) {
  //  for (int i = 0; i < m_items.size(); i++) m_items[i]->SetSelectable(e);
}

void Layer::selectAll() {
  for (int i = 0; i < m_items.size(); i++)
    if (m_items[i]) m_items[i]->setSelected(1);
}

void Layer::setBrushColor(QColor c) {
  m_brush.setColor(c);
  for (int i = 0; i < m_items.size(); i++)
    if (m_items[i]) m_items[i]->setBrush(m_brush);
}

void Layer::setPenColor(QColor c) {
  m_pen.setColor(c);
  for (int i = 0; i < m_items.size(); i++)
    if (m_items[i]) m_items[i]->setPen(m_pen);
}

Layer *Layer::loadFromXml(QXmlStreamReader *xml) {
  m_z = xml->attributes().value(tr("z")).toInt();
  m_pen.setColor(QColor(xml->attributes().value(tr("pen_color")).toUInt()));
  m_brush.setColor(QColor(xml->attributes().value(tr("brush_color")).toUInt()));
  xml->skipCurrentElement();
  return this;
}

bool Layer::saveToXml(QXmlStreamWriter *xml) {
  xml->writeStartElement(tr("layer"));
  xml->writeAttribute(tr("z"), QString("%1").arg(m_z));
  xml->writeAttribute(tr("pen_color"), QString("%1").arg(m_pen.color().rgba()));
  xml->writeAttribute(tr("brush_color"),
                      QString("%1").arg(m_brush.color().rgba()));
  xml->writeEndElement();
  return true;
}

GraphicsItem *Layer::selectItemByObjectName(QString name) {
  GraphicsItem *ret{nullptr};
  for (int i = 0; i < m_items.size(); i++) {
    if (m_items[i]->objectName() == name) {
      m_items[i]->setSelected(1);
      ret = m_items[i];
    } else
      m_items[i]->setSelected(0);
  }
  return ret;
}
