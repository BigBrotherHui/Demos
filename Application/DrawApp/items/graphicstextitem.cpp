#include "graphicstextitem.h"
#include <QGraphicsProxyWidget>
#include <QPlainTextEdit>
#include <QXmlStreamReader>
GraphicsTextItem::GraphicsTextItem(const QRect &rect, QGraphicsItem *parent)
    : GraphicsRectItem(rect, parent) {
  pMyProxy = new QGraphicsProxyWidget(this);
  QPlainTextEdit *ed = new QPlainTextEdit();
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
