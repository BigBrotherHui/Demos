#include "graphicsinstanceitem.h"
#include <QDebug>
#include <QGraphicsSceneEvent>
#include <QPainter>
#include "drawview.h"
#include "drawscene.h"
GraphicsInstanceItem::GraphicsInstanceItem(const QRect &rect,
                                           QGraphicsItem *parent)
    : GraphicsRectItem(rect, parent) {
      m_itemtype = instance;
    }

QRectF GraphicsInstanceItem::boundingRect() const { return rect(); }

QPainterPath GraphicsInstanceItem::shape() const {
    QPainterPath pa;
    pa.addRect(rect());
    return pa;
}

void GraphicsInstanceItem::paint(QPainter *painter,
                                 const QStyleOptionGraphicsItem *option,
                                 QWidget *widget) {
    if(pixmap.isNull())
    {
        painter->save();
        painter->fillRect(rect(),Qt::cyan);
        QFont ft("思源黑体 CN Medium",18);
        painter->setFont(ft);
        painter->setPen(Qt::red);
        painter->drawText(rect(),"instance",QTextOption(Qt::AlignCenter));
        painter->restore();
    }
    else
    {
        painter->drawPixmap(rect().toRect(), pixmap);
    }
}

QGraphicsItem *GraphicsInstanceItem::duplicate() const
{
    GraphicsInstanceItem *item = new GraphicsInstanceItem(rect().toRect(), dynamic_cast<QGraphicsItem *>(parent()));
    item->setPixmap(pixmap);
    item->createHandles();
    item->view=this->view;
    connect(view,&DrawView::updateScene,this,[=](QPixmap pix)
    {
        item->setPixmap(pix);
    });
    //    const_cast<GraphicsInstanceItem *>(this)->instanceItems.append(item);
    //    item->instanceItems.append(const_cast<GraphicsInstanceItem *>(this));
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
