#include "drawtool.h"
#include <QDebug>
#include <QGraphicsLineItem>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsView>
#include <QMenu>
#include <QMessageBox>
#include <QSet>
#include <QtMath>
#include <cmath>
#include "drawobj.h"
#include "globalsignalinstance.h"
#include "graphicsbezier.h"
#include "graphicslineitem.h"
#include "graphicstextitem.h"
#include "group.h"
#include "layer.h"
#include "layermanager.h"
#include "objectcontroller.h"
#include "mainwindow.h"
#define PI 3.1416

QList<DrawTool *> DrawTool::c_tools;
QPoint DrawTool::c_down;
QPoint DrawTool::c_last;
quint32 DrawTool::c_nDownFlags;

DrawShape DrawTool::c_drawShape = selection;

static SelectTool selectTool;
static RectTool rectTool(rectangle);
static RectTool roundRectTool(roundrect);
static RectTool ellipseTool(ellipse);
static RectTool textTool(text);
static RectTool instanceTool(instance);
static PolygonTool lineTool(line);
static PolygonTool polygonTool(polygon);
static PolygonTool bezierTool(bezier);
static PolygonTool polylineTool(polyline);
static RotationTool rotationTool;
enum SelectMode {
  none,
  netSelect,
  move,
  size,
  rotate,
  editor,
};

SelectMode selectMode = none;

int nDragHandle = Handle_None;

static void setCursor(DrawScene *scene, const QCursor &cursor) {
  QGraphicsView *view = scene->view();
  if (view) view->setCursor(cursor);
}

DrawTool::DrawTool(DrawShape shape) {
  m_drawShape = shape;
  m_hoverSizer = false;
  c_tools.push_back(this);
}

void DrawTool::mousePressEvent(QGraphicsSceneMouseEvent *event,
                               DrawScene *scene) {
    if(!MainWindow::currentView)
        return;
  if (!MainWindow::currentView->layerManager->isCurrentLayerValid()) {
//    QMessageBox::warning(nullptr, "warning", "please select layer first!");
    return;
  }
  auto p = event->scenePos().toPoint();
  QSize size = scene->gridSizeSpace();
  p.setX(p.x() / size.width() * size.width());
  p.setY(p.y() / size.height() * size.height());
  c_down = p;
  c_last = p;
  if (event->button() == Qt::RightButton && c_drawShape == selection) {
    if (!menu) {
      menu = new QMenu;
      QAction *addtogroup = new QAction("add to group");
      QAction *removegroup = new QAction("remove group");
      QAction *creategroup = new QAction("create a group with collidingItems");
      QAction *creategroupcascading =
          new QAction("create a group with cascading collidingItems");
      QAction *showallgroups = new QAction("show all groups");
      QObject::connect(creategroup, &QAction::triggered, this, [&] {
          QTransform t;
        GraphicsItem *selItem =
            dynamic_cast<GraphicsItem *>(scene->itemAt(event->scenePos(),t));
        if (!selItem) return;
        Group *g = new Group;
        QList<QGraphicsItem *> colItems = scene->collidingItems(selItem);
        for (int i = 0; i < colItems.size(); i++) {
          auto ii = dynamic_cast<GraphicsItem *>(colItems.at(i));
          if (ii) {
            g->addToGroup(ii);
          }
        }
        g->addToGroup(selItem);
        g->selectAll(1);
      });
      QObject::connect(creategroupcascading, &QAction::triggered, this, [&] {
          QTransform t;
        GraphicsItem *selItem =
            dynamic_cast<GraphicsItem *>(scene->itemAt(event->scenePos(),t));
        if (!selItem) return;
        Group *g = new Group;
        QSet<GraphicsItem *> items;
        scene->getAllCollidingItems(items, selItem);
        g->addToGroup(items);
        g->selectAll(1);
      });

      QObject::connect(showallgroups, &QAction::triggered, this, [&] {
          QTransform t;
        GraphicsItem *selItem =
            dynamic_cast<GraphicsItem *>(scene->itemAt(event->scenePos(),t));
        if (!selItem) return;
        for (auto g : selItem->groups()) {
          //          qDebug() << g->objectName() << g->items().size();
        }
        //        ObjectController::instance()->setObject(*selItem->groups().begin());
      });
      menu->addAction(creategroup);
      menu->addAction(creategroupcascading);
      menu->addAction(showallgroups);
    }
    menu->exec(event->screenPos());
  }
}

void DrawTool::mouseMoveEvent(QGraphicsSceneMouseEvent *event,
                              DrawScene *scene) {
  auto p = event->scenePos().toPoint();
  QSize size = scene->gridSizeSpace();
  p.setX(p.x() / size.width() * size.width());
  p.setY(p.y() / size.height() * size.height());
  c_last = p;
}

void DrawTool::mouseReleaseEvent(QGraphicsSceneMouseEvent *event,
                                 DrawScene *scene) {
  if (event->scenePos() == c_down)
    ;
  //    c_drawShape = selection;
  setCursor(scene, Qt::ArrowCursor);
}

void DrawTool::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event,
                                     DrawScene *scene) {}

DrawTool *DrawTool::findTool(DrawShape drawShape) {
  QList<DrawTool *>::const_iterator iter = c_tools.constBegin();
  for (; iter != c_tools.constEnd(); ++iter) {
    if ((*iter)->m_drawShape == drawShape) return (*iter);
  }
  return 0;
}

SelectTool::SelectTool() : DrawTool(selection) {
  dashRect = 0;
  selLayer = 0;
  opposite_ = QPointF();
}

void SelectTool::mousePressEvent(QGraphicsSceneMouseEvent *event,
                                 DrawScene *scene) {
  DrawTool::mousePressEvent(event, scene);

  if (event->button() != Qt::LeftButton) return;

  if (!m_hoverSizer) scene->mouseEvent(event);

  nDragHandle = Handle_None;
  selectMode = none;
  QList<QGraphicsItem *> items = scene->selectedItems();
  AbstractShape *item = 0;

  if (items.count() == 1)
    item = qgraphicsitem_cast<AbstractShape *>(items.first());

  if (item != 0) {
    nDragHandle = item->collidesWithHandle(event->scenePos());
    if (nDragHandle != Handle_None && nDragHandle <= Left)
      selectMode = size;
    else if (nDragHandle > Left)
      selectMode = editor;
    else
      selectMode = move;

    if (nDragHandle != Handle_None && nDragHandle <= Left) {
      opposite_ = item->opposite(nDragHandle);
      if (opposite_.x() == 0) opposite_.setX(1);
      if (opposite_.y() == 0) opposite_.setY(1);
    }

    setCursor(scene, Qt::ClosedHandCursor);

  } else if (items.count() > 1)
    selectMode = move;

  if (selectMode == none) {
    selectMode = netSelect;
    if (scene->view()) {
      QGraphicsView *view = scene->view();
      view->setDragMode(QGraphicsView::RubberBandDrag);
    }
#if 0
        if ( selLayer ){
            scene->destroyGroup(selLayer);
            selLayer = 0;
        }
#endif
  }

  if (selectMode == move && items.count() == 1) {
    if (dashRect) {
      scene->removeItem(dashRect);
      delete dashRect;
      dashRect = 0;
    }

    dashRect = new QGraphicsPathItem(item->shape());
    dashRect->setPen(Qt::DashLine);
    dashRect->setPos(item->pos());
    dashRect->setTransformOriginPoint(item->transformOriginPoint());
    dashRect->setTransform(item->transform());
    dashRect->setRotation(item->rotation());
    dashRect->setScale(item->scale());
    dashRect->setZValue(item->zValue());
    scene->addItem(dashRect);
    initialPositions = item->pos();
  }
}

void SelectTool::mouseMoveEvent(QGraphicsSceneMouseEvent *event,
                                DrawScene *scene) {
  DrawTool::mouseMoveEvent(event, scene);
  QList<QGraphicsItem *> items = scene->selectedItems();
  AbstractShape *item = 0;
  if (items.count() == 1) {
    item = qgraphicsitem_cast<AbstractShape *>(items.first());
    if (item != 0) {
      if (nDragHandle != Handle_None && selectMode == size) {
        if (opposite_.isNull()) {
          opposite_ = item->opposite(nDragHandle);
          if (opposite_.x() == 0) opposite_.setX(1);
          if (opposite_.y() == 0) opposite_.setY(1);
        }

        QPointF new_delta = item->mapFromScene(c_last) - opposite_;
        QPointF initial_delta = item->mapFromScene(c_down) - opposite_;

        double sx = new_delta.x() / initial_delta.x();
        double sy = new_delta.y() / initial_delta.y();

        item->stretch(nDragHandle, sx, sy, opposite_);

        emit scene->itemResize(item, nDragHandle, QPointF(sx, sy));

        //  qDebug()<<"scale:"<<nDragHandle<< item->mapToScene(opposite_)<< sx
        //  << " ，" << sy
        //         << new_delta << item->mapFromScene(c_last)
        //         << initial_delta << item->mapFromScene(c_down) <<
        //         item->boundingRect();

      } else if (nDragHandle > Left && selectMode == editor) {
        item->control(nDragHandle, c_last);
        emit scene->itemControl(item, nDragHandle, c_last, c_down);
      } else if (nDragHandle == Handle_None) {
        auto p = event->scenePos().toPoint();
        QSize size = scene->gridSizeSpace();
        p.setX(p.x() / size.width() * size.width());
        p.setY(p.y() / size.height() * size.height());
        int handle = item->collidesWithHandle(p);
        if (handle != Handle_None) {
          setCursor(scene, Qt::OpenHandCursor);
          m_hoverSizer = true;
        } else {
          setCursor(scene, Qt::ArrowCursor);
          m_hoverSizer = false;
        }
      }
    }
  }

  if (selectMode == move) {
    setCursor(scene, Qt::ClosedHandCursor);
    if (dashRect) {
        dashRect->setPos(initialPositions + c_last - c_down);
    }
  }

  if (selectMode != size && items.count() > 1) {
    scene->mouseEvent(event);
  }
}

void SelectTool::mouseReleaseEvent(QGraphicsSceneMouseEvent *event,
                                   DrawScene *scene) {
  DrawTool::mouseReleaseEvent(event, scene);

  if (event->button() != Qt::LeftButton) return;

  QList<QGraphicsItem *> items = scene->selectedItems();
  if (items.count() == 1) {
    AbstractShape *item = qgraphicsitem_cast<AbstractShape *>(items.first());
    if (item != 0 && selectMode == move && c_last != c_down) {
      item->setPos(initialPositions + c_last - c_down);
      emit scene->itemMoved(item, c_last - c_down);
    } else if (item != 0 && (selectMode == size || selectMode == editor) &&
               c_last != c_down) {
      item->updateCoordinate();
    }
  } else if (items.count() > 1 && selectMode == move && c_last != c_down) {
    emit scene->itemMoved(NULL, c_last - c_down);
  }

  if (selectMode == netSelect) {
    if (scene->view()) {
      QGraphicsView *view = scene->view();
      view->setDragMode(QGraphicsView::NoDrag);
    }
#if 0
        if ( scene->selectedItems().count() > 1 ){
            selLayer = scene->createGroup(scene->selectedItems());
            selLayer->setSelected(true);
        }
#endif
  }
  if (dashRect) {
    scene->removeItem(dashRect);
    delete dashRect;
    dashRect = 0;
  }
  selectMode = none;
  nDragHandle = Handle_None;
  m_hoverSizer = false;
  opposite_ = QPointF();
  scene->mouseEvent(event);
}

RotationTool::RotationTool() : DrawTool(rotation) {
  lastAngle = 0;
  dashRect = 0;
}

void RotationTool::mousePressEvent(QGraphicsSceneMouseEvent *event,
                                   DrawScene *scene) {
  DrawTool::mousePressEvent(event, scene);
  if (event->button() != Qt::LeftButton) return;

  if (!m_hoverSizer) scene->mouseEvent(event);

  QList<QGraphicsItem *> items = scene->selectedItems();
  if (items.count() == 1) {
    AbstractShape *item = qgraphicsitem_cast<AbstractShape *>(items.first());
    if (item != 0) {
      nDragHandle = item->collidesWithHandle(event->scenePos());
      if (nDragHandle != Handle_None) {
        QPointF origin = item->mapToScene(item->boundingRect().center());

        qreal len_y = c_last.y() - origin.y();
        qreal len_x = c_last.x() - origin.x();

        qreal angle = atan2(len_y, len_x) * 180 / PI;

        lastAngle = angle;
        selectMode = rotate;

        if (dashRect) {
          scene->removeItem(dashRect);
          delete dashRect;
          dashRect = 0;
        }

        dashRect = new QGraphicsPathItem(item->shape());
        dashRect->setPen(Qt::DashLine);
        dashRect->setPos(item->pos());
        dashRect->setTransformOriginPoint(item->transformOriginPoint());
        dashRect->setTransform(item->transform());
        dashRect->setRotation(item->rotation());
        dashRect->setScale(item->scale());
        dashRect->setZValue(item->zValue());
        scene->addItem(dashRect);
        setCursor(scene, QCursor((QPixmap(":/icons/rotate.png"))));
      } else {
        c_drawShape = selection;
        selectTool.mousePressEvent(event, scene);
      }
    }
  }
}

void RotationTool::mouseMoveEvent(QGraphicsSceneMouseEvent *event,
                                  DrawScene *scene) {
  DrawTool::mouseMoveEvent(event, scene);
  QList<QGraphicsItem *> items = scene->selectedItems();
  if (items.count() == 1) {
    AbstractShape *item = qgraphicsitem_cast<AbstractShape *>(items.first());
    if (item != 0 && nDragHandle != Handle_None && selectMode == rotate) {
      QPointF origin = item->mapToScene(item->boundingRect().center());

      qreal len_y = c_last.y() - origin.y();
      qreal len_x = c_last.x() - origin.x();
      qreal angle = atan2(len_y, len_x) * 180 / PI;

      angle = item->rotation() + int(angle - lastAngle);

      if (angle > 360) angle -= 360;
      if (angle < -360) angle += 360;

      if (dashRect) {
        // dashRect->setTransform(QTransform::fromTranslate(15,15),true);
        // dashRect->setTransform(QTransform().rotate(angle));
        // dashRect->setTransform(QTransform::fromTranslate(-15,-15),true);
        dashRect->setRotation(angle);
      }

      setCursor(scene, QCursor((QPixmap(":/icons/rotate.png"))));
    } else if (item) {
      int handle = item->collidesWithHandle(event->scenePos());
      if (handle != Handle_None) {
        setCursor(scene, QCursor((QPixmap(":/icons/rotate.png"))));
        m_hoverSizer = true;
      } else {
        setCursor(scene, Qt::ArrowCursor);
        m_hoverSizer = false;
      }
    }
  }
  scene->mouseEvent(event);
}

void RotationTool::mouseReleaseEvent(QGraphicsSceneMouseEvent *event,
                                     DrawScene *scene) {
  DrawTool::mouseReleaseEvent(event, scene);
  if (event->button() != Qt::LeftButton) return;

  QList<QGraphicsItem *> items = scene->selectedItems();
  if (items.count() == 1) {
    AbstractShape *item = qgraphicsitem_cast<AbstractShape *>(items.first());
    if (item != 0 && nDragHandle != Handle_None && selectMode == rotate) {
      QPointF origin = item->mapToScene(item->boundingRect().center());
      QPointF delta = c_last - origin;
      qreal len_y = c_last.y() - origin.y();
      qreal len_x = c_last.x() - origin.x();
      qreal angle = atan2(len_y, len_x) * 180 / PI, oldAngle = item->rotation();
      angle = item->rotation() + int(angle - lastAngle);

      if (angle > 360) angle -= 360;
      if (angle < -360) angle += 360;

      item->setRotation(angle);
      emit scene->itemRotate(item, oldAngle);
      qDebug() << "rotate:" << angle << item->boundingRect();
    }
  }

  setCursor(scene, Qt::ArrowCursor);
  selectMode = none;
  nDragHandle = Handle_None;
  lastAngle = 0;
  m_hoverSizer = false;
  if (dashRect) {
    scene->removeItem(dashRect);
    delete dashRect;
    dashRect = 0;
  }
  scene->mouseEvent(event);
}

RectTool::RectTool(DrawShape drawShape) : DrawTool(drawShape) {
  item = 0;
  connect(this, &RectTool::signal_instanceItemDoubleClicked,
          GlobalSignalInstance::instance(),
          &GlobalSignalInstance::signal_wrapper_instanceItemDoubleClicked);
}

void RectTool::mousePressEvent(QGraphicsSceneMouseEvent *event,
                               DrawScene *scene) {
  if (event->button() != Qt::LeftButton) return;
  scene->clearSelection();
  DrawTool::mousePressEvent(event, scene);
  switch (c_drawShape) {
    case rectangle:
      item = new GraphicsRectItem(QRect(1, 1, 1, 1));
      break;
    case roundrect:
      item = new GraphicsRectItem(QRect(1, 1, 1, 1), true);
      break;
      //    case ellipse:
      //      item = new GraphicsEllipseItem(QRect(1, 1, 1, 1));
      //      break;
    case text:
      item = new GraphicsTextItem(QRect(1, 1, 1, 1));
      break;
    case instance:
      item = new GraphicsInstanceItem(QRect(1, 1, 1, 1));
  }
  if (item == 0) return;
  item->createHandles();
  c_down += QPoint(2, 2);
  auto p = event->scenePos().toPoint();
  QSize sz = scene->gridSizeSpace();
  p.setX(p.x() / sz.width() * sz.width());
  p.setY(p.y() / sz.height() * sz.height());
  item->setPos(p);
  if(MainWindow::currentView)
      MainWindow::currentView->layerManager->GetCurrentLayer()->AddItem(item);
  item->setSelected(true);

  selectMode = size;
  nDragHandle = RightBottom;
}

void RectTool::mouseMoveEvent(QGraphicsSceneMouseEvent *event,
                              DrawScene *scene) {
  setCursor(scene, Qt::CrossCursor);
  selectTool.mouseMoveEvent(event, scene);
}

void RectTool::mouseReleaseEvent(QGraphicsSceneMouseEvent *event,
                                 DrawScene *scene) {
  selectTool.mouseReleaseEvent(event, scene);
  auto p = event->scenePos().toPoint();
  QSize size = scene->gridSizeSpace();
  p.setX(p.x() / size.width() * size.width());
  p.setY(p.y() / size.height() * size.height());
  if (p == (c_down - QPoint(2, 2))) {
    if (item != 0) {
      item->setSelected(false);
      if(MainWindow::currentView)
          MainWindow::currentView->layerManager->GetCurrentLayer()->RemoveItem(item);
      //      emit scene->itemRemoved();
      item = 0;
    }
    qDebug() << "RectTool removeItem:";
  } else if (item) {
    emit scene->itemAdded(item);
    QPoint r = item->pos().toPoint();
    static_cast<GraphicsRectItem *>(item)->setPos(
        QPoint(r.x() / size.width() * size.width(),
               r.y() / size.height() * size.height()));
    item->setWidth(qCeil(item->width() + size.width() / 2) / size.width() *
                   size.width());
    item->setHeight(qCeil(item->height() + size.width() / 2) / size.height() *
                    size.height());
    if (dynamic_cast<GraphicsTextItem *>(item))
      c_drawShape = selection;
    else if (dynamic_cast<GraphicsInstanceItem *>(item))
      c_drawShape = selection;
  }
}

void RectTool::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event,
                                     DrawScene *scene) {
  if (event->button() == Qt::LeftButton) {
    if (!item) {
      item = new GraphicsInstanceItem(QRect(1, 1, 1, 1));
      if(MainWindow::currentView)
          MainWindow::currentView->layerManager->GetCurrentLayer()->AddItem(item);
    }
    connect(this, &RectTool::signal_instanceItemDoubleClicked,
            GlobalSignalInstance::instance(),
            &GlobalSignalInstance::signal_wrapper_instanceItemDoubleClicked,
            static_cast<Qt::ConnectionType>(Qt::DirectConnection |
                                            Qt::UniqueConnection));
    emit signal_instanceItemDoubleClicked(
        static_cast<GraphicsInstanceItem *>(item));
  }
}

PolygonTool::PolygonTool(DrawShape shape) : DrawTool(shape) {
  item = NULL;
  m_nPoints = 0;
}

void PolygonTool::mousePressEvent(QGraphicsSceneMouseEvent *event,
                                  DrawScene *scene) {
  DrawTool::mousePressEvent(event, scene);
  auto p = event->scenePos().toPoint();
  QSize sz = scene->gridSizeSpace();
  p.setX(p.x() / sz.width() * sz.width());
  p.setY(p.y() / sz.height() * sz.height());
  if (event->button() == Qt::RightButton) {
    if (!item) return;
    item->endPoint(p);
    item->updateCoordinate();
    emit scene->itemAdded(item);
    item = NULL;
    selectMode = none;
    //    c_drawShape = selection;
    m_nPoints = 0;
  }
  if (event->button() != Qt::LeftButton) return;

  if (item == NULL) {
    if (c_drawShape == polygon) {
      item = new GraphicsPolygonItem(NULL);
    } else if (c_drawShape == bezier) {
      item = new GraphicsBezier();
    } else if (c_drawShape == polyline) {
      item = new GraphicsBezier(false);
    } else if (c_drawShape == line) {
      item = new GraphicsLineItem(0);
    }
    item->createHandles();
    item->setPos(p);
    if(MainWindow::currentView)
        MainWindow::currentView->layerManager->GetCurrentLayer()->AddItem(item);
    initialPositions = c_down;
    item->addPoint(c_down);
    item->setSelected(true);
    m_nPoints++;

  } else if (c_down == c_last) {
    /*
    if ( item != NULL )
    {
        scene->removeItem(item);
        delete item;
        item = NULL ;
        c_drawShape = selection;
        selectMode = none;
        return ;
    }
    */
  }
  item->addPoint(c_down + QPoint(1, 0));
  m_nPoints++;
  selectMode = size;
  nDragHandle = item->handleCount();
}

void PolygonTool::mouseMoveEvent(QGraphicsSceneMouseEvent *event,
                                 DrawScene *scene) {
  DrawTool::mouseMoveEvent(event, scene);
  setCursor(scene, Qt::CrossCursor);

  //    selectTool.mouseMoveEvent(event,scene);

  if (item != 0) {
    if (nDragHandle != Handle_None && selectMode == size) {
      item->control(nDragHandle, c_last);
    }
  }
}

void PolygonTool::mouseReleaseEvent(QGraphicsSceneMouseEvent *event,
                                    DrawScene *scene) {
  DrawTool::mousePressEvent(event, scene);
  if (c_drawShape == line) {
    auto p = event->scenePos().toPoint();
    QSize size = scene->gridSizeSpace();
    p.setX(p.x() / size.width() * size.width());
    p.setY(p.y() / size.height() * size.height());
    item->endPoint(p);
    item->updateCoordinate();
    emit scene->itemAdded(item);
    item = NULL;
    selectMode = none;
    //    c_drawShape = selection;
    m_nPoints = 0;
  }
}

void PolygonTool::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event,
                                        DrawScene *scene) {
  DrawTool::mouseDoubleClickEvent(event, scene);
  //  item->endPoint(event->scenePos());
  //  item->updateCoordinate();
  //  emit scene->itemAdded(item);
  //  item = NULL;
  //  selectMode = none;
  //  c_drawShape = selection;
  //  m_nPoints = 0;
}
