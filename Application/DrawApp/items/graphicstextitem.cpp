#include "graphicstextitem.h"
#include <QGraphicsProxyWidget>
#include <QPlainTextEdit>
#include <QXmlStreamReader>
#include "drawscene.h"
GraphicsTextItem::GraphicsTextItem(const QRect &rect, QGraphicsItem *parent)
    : GraphicsRectItem(rect, parent) {
  pMyProxy = new QGraphicsProxyWidget(this);
  m_itemtype = text;
  ed = new QPlainTextEdit();
  ed->setMinimumSize(QSize(1, 1));
  ed->setPlaceholderText("input text");
  pMyProxy->setWidget(ed);
  pMyProxy->setPos(boundingRect().center() - ed->rect().center());
}

void GraphicsTextItem::paint(QPainter *painter,
                             const QStyleOptionGraphicsItem *option,
                             QWidget *widget) {
  pMyProxy->setGeometry(rect().adjusted(1, 1, -1, -1));
}

bool GraphicsTextItem::loadFromXml(QXmlStreamReader *xml) {
  readBaseAttributes(xml);
  updateCoordinate();
  static_cast<QPlainTextEdit *>(pMyProxy->widget())
      ->setPlainText(xml->attributes().value(tr("text")).toString());
  pMyProxy->setGeometry(rect().adjusted(1, 1, -1, -1));
  xml->skipCurrentElement();
  return 1;
}

bool GraphicsTextItem::saveToXml(QXmlStreamWriter *xml) {
  xml->writeStartElement(tr("text"));
  writeBaseAttributes(xml);
  xml->writeAttribute(
      "text", static_cast<QPlainTextEdit *>(pMyProxy->widget())->toPlainText());
  xml->writeEndElement();
  return 1;
}

QGraphicsItem *GraphicsTextItem::duplicate() const
{
    GraphicsTextItem *item = new GraphicsTextItem(rect().toRect(), dynamic_cast<QGraphicsItem *>(parent()));
    item->ed->setPlainText(ed->toPlainText());
    item->createHandles();
    item->m_width = width();
    item->m_height = height();
    item->setPos(pos().x(), pos().y());
    item->setPen(pen());
    item->setBrush(brush());
    item->setTransform(transform());
    item->setTransformOriginPoint(transformOriginPoint());
    item->setRotation(rotation());
    item->setScale(scale());
    item->setZValue(zValue());
    item->m_fRatioY = m_fRatioY;
    item->m_fRatioX = m_fRatioX;
    item->updateCoordinate();
    return item;
}
