#include "drawview.h"
#include <QProcess>
#include <QSvgGenerator>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>
#include <QtConcurrent/QtConcurrent>
#include "drawscene.h"
#include "gdsimporter.h"
#include "graphicsbezier.h"
#include "graphicslineitem.h"
#include "graphicsrectitem.h"
#include "graphicstextitem.h"
#include "layer.h"
#include "layermanager.h"
#include "ruler.h"
#include "mainwindow.h"
#include "mapview.h"
#include "gdsimporter.h"
DrawView::DrawView(QGraphicsScene *pscene) : QGraphicsView(pscene) {
    layerManager=new LayerManager(dynamic_cast<DrawScene *>(pscene),this);
  m_hruler = new QtRuleBar(Qt::Horizontal, this, this);
  m_vruler = new QtRuleBar(Qt::Vertical, this, this);
  box = new QtCornerBox(this);
  setRenderHint(QPainter::Antialiasing, false);
  //  setCacheMode(QGraphicsView::CacheBackground);
  setOptimizationFlags(QGraphicsView::DontSavePainterState);
  setViewportUpdateMode(QGraphicsView::SmartViewportUpdate);
  // view->setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
  setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
  setViewport(new QWidget);

//  setAttribute(Qt::WA_DeleteOnClose);
  isUntitled = true;

  modified = false;

  //  this->setMouseTracking(true);

  this->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  this->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  //  this->verticalScrollBar()->blockSignals(true);
  //  this->horizontalScrollBar()->blockSignals(true);

  //  this->horizontalScrollBar()->setEnabled(false);
  //  this->verticalScrollBar()->setEnabled(false);
  //  this->verticalScrollBar()->setValue(this->verticalScrollBar()->minimum());
  //  this->horizontalScrollBar()->setValue(this->horizontalScrollBar()->minimum());
  ruler = new Ruler(this);
  connect(ruler, &Ruler::signal_zoomFactorChanged, this, [&](int v) {
    qreal scale = qPow(qreal(2), (v - 250) / qreal(50));
    QMatrix matrix;
    matrix.scale(scale, scale);
    setMatrix(matrix);
    updateRuler();
  });

  mapView= new MapView(this, this);
  connect(static_cast<DrawScene *>(scene()), &DrawScene::itemAdded, mapView, &MapView::updateImage);
  connect(static_cast<DrawScene *>(scene()), &DrawScene::itemRemoved, mapView, &MapView::updateImage);
  connect(static_cast<DrawScene *>(scene()), &DrawScene::itemResize, mapView, &MapView::updateImage);
  connect(static_cast<DrawScene *>(scene()), &DrawScene::itemMoved, mapView, &MapView::updateImage);
  connect(this, &DrawView::sizeChanged, mapView, &MapView::updateImage);
}

DrawView::~DrawView()
{
    qDebug()<<"~DrawView";
}

void DrawView::zoomIn() { ruler->updateValue(6); }

void DrawView::zoomOut() { ruler->updateValue(-6); }

void DrawView::newFile() {
  static int sequenceNumber = 1;

  isUntitled = true;
  curFile = tr("drawboard%1.gds").arg(sequenceNumber++);
  setWindowTitle(curFile + "[*]");
}

bool DrawView::loadFile(const QString &fileName) {
  QFile file(fileName);
  if (!file.open(QFile::ReadOnly | QFile::Text)) {
    QMessageBox::warning(
        this, tr("EDA DrawBoard"),
        tr("Cannot read file %1:\n%2.").arg(fileName).arg(file.errorString()));
    return false;
  }
  scene()->setSceneRect(0, 0, 30000, 30000);
  QProcess process;
  QStringList arguments;
  arguments << "-i" << fileName;
  process.setProgram("gdscpp.exe");
  process.setArguments(arguments);
  process.setStandardOutputFile("output.txt");

  process.start();
  process.waitForFinished();
  loadCanvas("output.txt");
  setCurrentFile(fileName);
  QFile::remove("output.txt");
  return 1;
}

bool DrawView::save()
{
    if (isUntitled) {
        return saveAs();
    } else {
        return saveFile(curFile);
    }
}

bool DrawView::saveAs() {
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save As"), curFile);
  if (fileName.isEmpty()) return false;

  return saveFile(fileName);
}

bool DrawView::saveFile(const QString &fileName) {
  QFile file(fileName);
  if (!file.open(QFile::WriteOnly | QFile::Text)) {
    QMessageBox::warning(
        this, tr("Qt Drawing"),
        tr("Cannot write file %1:\n%2.").arg(fileName).arg(file.errorString()));
    return false;
  }
  std::vector<GraphicsItem *> items;
  auto layers=layerManager->GetAllLayers();
  for(int i=0;i<layers.size();i++)
  {
      auto its=layers[i]->GetAllItems();
      for(int k=0;k<its.size();k++)
          items.push_back(its[k]);
  }
  gdsimporter::instance()->writeFile(fileName,items);
//  QXmlStreamWriter xml(&file);
//  xml.setAutoFormatting(true);
//  xml.writeStartDocument();
//  xml.writeDTD("<!DOCTYPE DrawApp>");
//  xml.writeStartElement("canvas");
//  xml.writeAttribute("width", QString("%1").arg(scene()->width()));
//  xml.writeAttribute("height", QString("%1").arg(scene()->height()));
//  foreach (QGraphicsItem *item, scene()->items()) {
//    AbstractShape *ab = qgraphicsitem_cast<AbstractShape *>(item);
//    QGraphicsItemGroup *g =
//        dynamic_cast<QGraphicsItemGroup *>(item->parentItem());
//    if (ab && !qgraphicsitem_cast<SizeHandleRect *>(ab) && !g) {
//      ab->saveToXml(&xml);
//    }
//  }
//  if(MainWindow::currentView)
//  {
//      auto alllayers=MainWindow::currentView->layerManager->GetAllLayers();
//      for (Layer *l : alllayers) {
//          l->saveToXml(&xml);
//      }
//  }
//  xml.writeEndElement();
//  xml.writeEndDocument();
#if 0
    QSvgGenerator generator;
    generator.setFileName(fileName);
    generator.setSize(QSize(800, 600));
    generator.setTitle(tr("SVG Generator Example Drawing"));
    generator.setDescription(tr("An SVG drawing created by the SVG Generator "
                                "Example provided with Qt."));
//![configure SVG generator]
//![begin painting]
    QPainter painter;
    painter.begin(&generator);
//![begin painting]
//!
    scene()->clearSelection();
    scene()->render(&painter);
//![end painting]
    painter.end();
//![end painting]
#endif
  setCurrentFile(fileName);
  return true;
}

QString DrawView::userFriendlyCurrentFile() { return strippedName(curFile); }

void DrawView::show()
{
    m_isSaved=0;
    return QGraphicsView::show();
}

void DrawView::closeEvent(QCloseEvent *event) {
  if (!m_isSaved) {
      if(maybeSave())
      {
          event->accept();
      }
      else {
          event->ignore();
      }
      m_isSaved=1;
  }
  else
      event->accept();
}

void DrawView::wheelEvent(QWheelEvent *event) {
  //  if (event->modifiers() == Qt::CTRL) {
  if ((event->delta() > 0) && (m_scale >= 50)) {
    return;
  } else if ((event->delta() < 0) && (m_scale <= 0.01)) {
    return;
  } else {
    qreal scaleFactor = this->matrix().m11();
    m_scale = scaleFactor;

    int wheelDeltaValue = event->delta();

    if (wheelDeltaValue > 0) {
      //      this->scale(1.2, 1.2);
      ruler->updateValue(6);
    } else {
      //      this->scale(1.0 / 1.2, 1.0 / 1.2);
      ruler->updateValue(-6);
    }
    //    updateRuler();
  }
  //}
}

void DrawView::mouseMoveEvent(QMouseEvent *event) {
  QPointF pt = mapToScene(event->pos());
  m_hruler->updatePosition(event->pos());
  m_vruler->updatePosition(event->pos());
  emit positionChanged(pt.x(), pt.y());
  QGraphicsView::mouseMoveEvent(event);
}

void DrawView::resizeEvent(QResizeEvent *event) {
  QGraphicsView::resizeEvent(event);
  this->setViewportMargins(RULER_SIZE - 1, RULER_SIZE - 1, 0, 0);
  m_hruler->resize(this->size().width() - RULER_SIZE - 1, RULER_SIZE);
  m_hruler->move(RULER_SIZE, 0);
  m_vruler->resize(RULER_SIZE, this->size().height() - RULER_SIZE - 1);
  m_vruler->move(0, RULER_SIZE);

  box->resize(RULER_SIZE, RULER_SIZE);
  box->move(0, 0);
  updateRuler();
  ruler->move(20, height() - ruler->height() - 20);
}

void DrawView::scrollContentsBy(int dx, int dy) {
  QGraphicsView::scrollContentsBy(dx, dy);
  updateRuler();
}

void DrawView::updateRuler() {
  if (scene() == 0) return;
  QRectF viewbox = this->rect();
  QPointF offset = mapFromScene(scene()->sceneRect().topLeft());
  double factor = 1. / transform().m11();
  double lower_x = factor * (viewbox.left() - offset.x());
  double upper_x = factor * (viewbox.right() - RULER_SIZE - offset.x());
  m_hruler->setRange(lower_x, upper_x, upper_x - lower_x);
  m_hruler->update();

  double lower_y = factor * (viewbox.top() - offset.y()) * -1;
  double upper_y = factor * (viewbox.bottom() - RULER_SIZE - offset.y()) * -1;

  m_vruler->setRange(lower_y, upper_y, upper_y - lower_y);
  m_vruler->update();
  emit sizeChanged();
  // qDebug()<<viewbox<<QPoint(lower_x,upper_x) << QPoint(lower_y,upper_y) <<
  // offset;
}

bool DrawView::maybeSave() {
  if (isModified()) {
    QMessageBox::StandardButton ret;
    ret = QMessageBox::warning(
        this, tr("MDI"),
        tr("'%1' has been modified.\n"
           "Do you want to save your changes?")
            .arg(userFriendlyCurrentFile()),
        QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
    if (ret == QMessageBox::Save)
    {
        emit updateScene(grab());
        if(m_isSubView)
        {
            return 1;
        }
        return save();
    }
    else if (ret == QMessageBox::Cancel)
      return false;
  }
  return true;
}

void DrawView::setCurrentFile(const QString &fileName) {
  curFile = QFileInfo(fileName).canonicalFilePath();
  isUntitled = false;
  setModified(false);
  setWindowModified(false);
  setWindowTitle(userFriendlyCurrentFile() + "[*]");
}

QString DrawView::strippedName(const QString &fullFileName) {
  return QFileInfo(fullFileName).fileName();
}

void DrawView::loadCanvas(QXmlStreamReader *xml) {
  Q_ASSERT(xml->isStartElement() && xml->name() == "canvas");
  std::vector<GraphicsItem *> items;
  if(!MainWindow::currentView)
      return;
  while (xml->readNextStartElement()) {
    if (xml->name() == tr("layer")) {
      Layer *l = MainWindow::currentView->layerManager->loadFromXml(
          xml, dynamic_cast<DrawScene *>(scene()));
      MainWindow::currentView->layerManager->AddLayer(l);
    }
    GraphicsItem *item = NULL;
    if (xml->name() == tr("rect")) {
      item = new GraphicsRectItem(QRect(0, 0, 1, 1));
    } else if (xml->name() == tr("roundrect")) {
      item = new GraphicsRectItem(QRect(0, 0, 1, 1), true);
    } /*else if (xml->name() == tr("ellipse"))
      item = new GraphicsEllipseItem(QRect(0, 0, 1, 1));*/
    else if (xml->name() == tr("polygon"))
      item = new GraphicsPolygonItem();
    else if (xml->name() == tr("bezier"))
      item = new GraphicsBezier();
    else if (xml->name() == tr("polyline"))
      item = new GraphicsBezier(false);
    else if (xml->name() == tr("line"))
      item = new GraphicsLineItem();
    else if (xml->name() == "text")
      item = new GraphicsTextItem(QRect(0, 0, 1, 1));
    //    else if (xml->name() == tr("group"))
    //      item = qgraphicsitem_cast<AbstractShape *>(loadGroupFromXML(xml));
    //    else
    //      xml->skipCurrentElement();
    item->createHandles();
    if (item && item->loadFromXml(xml)) {
      items.push_back(item);
    } else if (item)
      delete item;
  }
  for (int i = 0; i < items.size(); i++) {
    auto item = items[i];
    MainWindow::currentView->layerManager->AddToSpecifiedLayer(item->zValue(), item);
    scene()->addItem(item);
  }
}

void DrawView::loadCanvas(QString gdsfilepath) {
  //  LayerManager::GetInstance()->AddLayer(dynamic_cast<DrawScene *>(scene()));
  std::vector<GraphicsItem *> items;
  QFuture<void> f = QtConcurrent::run(
      [&] { gdsimporter::instance()->parseFile(gdsfilepath, items); });
  QFutureWatcher<void> w;
  w.setFuture(f);
  while (w.isRunning()) {
    qApp->processEvents(QEventLoop::AllEvents);
  }
  for (int i = 0; i < items.size(); i++) {
    auto item = items[i];
    //    LayerManager::GetInstance()->AddToSpecifiedLayer(0, item);
    scene()->addItem(item);
  }
}

GraphicsItemGroup *DrawView::loadGroupFromXML(QXmlStreamReader *xml) {
  QList<QGraphicsItem *> items;
  qreal angle = xml->attributes().value(tr("rotate")).toDouble();
  while (xml->readNextStartElement()) {
    AbstractShape *item = NULL;
    if (xml->name() == tr("rect")) {
      item = new GraphicsRectItem(QRect(0, 0, 1, 1));
    } else if (xml->name() == tr("roundrect")) {
      item = new GraphicsRectItem(QRect(0, 0, 1, 1), true);
    } /*else if (xml->name() == tr("ellipse"))
      item = new GraphicsEllipseItem(QRect(0, 0, 1, 1));*/
    else if (xml->name() == tr("polygon"))
      item = new GraphicsPolygonItem();
    else if (xml->name() == tr("bezier"))
      item = new GraphicsBezier();
    else if (xml->name() == tr("polyline"))
      item = new GraphicsBezier(false);
    else if (xml->name() == tr("line"))
      item = new GraphicsLineItem();
    else if (xml->name() == tr("group"))
      item = qgraphicsitem_cast<AbstractShape *>(loadGroupFromXML(xml));
    else if (xml->name() == "text")
      item = new GraphicsTextItem(QRect(0, 0, 1, 1));
    else
      xml->skipCurrentElement();
    if (item && item->loadFromXml(xml)) {
      scene()->addItem(item);
      items.append(item);
    } else if (item)
      delete item;
  }

  if (items.count() > 0) {
    DrawScene *s = dynamic_cast<DrawScene *>(scene());
    GraphicsItemGroup *group = s->createGroup(items, false);
    if (group) {
      group->setRotation(angle);
      group->updateCoordinate();
      // qDebug()<<"angle:" <<angle;
    }
    return group;
  }
  return 0;
}
