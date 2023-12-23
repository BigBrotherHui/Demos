#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QMainWindow>
#include <QMap>
#include <QTimer>
#include "drawscene.h"
#include "drawview.h"
#include "objectcontroller.h"

QT_BEGIN_NAMESPACE
class QAction;
class QActionGroup;
class QToolBox;
class QSpinBox;
class QComboBox;
class QFontComboBox;
class QButtonGroup;
class QLineEdit;
class QGraphicsTextItem;
class QFont;
class QToolButton;
class QAbstractButton;
class QGraphicsView;
class QListWidget;
class QStatusBar;
class QMdiArea;
class QMdiSubWindow;
class QSignalMapper;
class QUndoStack;
class QUndoView;
class CommandWidget;
QT_END_NAMESPACE

class QtVariantProperty;
class QtProperty;
class Group;
class MainWindow : public QMainWindow {
  Q_OBJECT

 public:
  MainWindow(QWidget *parent = 0);
  ~MainWindow();

  bool openFile(const QString &fileName);

 public slots:

  void newFile();
  void open();
  void save();
  DrawView *createMdiChild();
  void updateMenus();
  void updateWindowMenu();
  void setActiveSubWindow(QWidget *window);

  void addShape();
  void updateActions();
  void deleteItem();
  void itemSelected();
  void itemMoved(QGraphicsItem *item, const QPointF &oldPosition);
  void itemAdded(QGraphicsItem *item);
  void itemRotate(QGraphicsItem *item, const qreal oldAngle);
  void itemResize(QGraphicsItem *item, int handle, const QPointF &scale);
  void itemControl(QGraphicsItem *item, int handle, const QPointF &newPos,
                   const QPointF &lastPos_);
  void groupChanged(Group *g);
  void addLayer();
  void removeLayer();
  void layerChanged();
  void on_actionBringToFront_triggered();
  void on_actionSendToBack_triggered();
  void on_aglin_triggered();
  void zoomIn();
  void zoomOut();
  void on_group_triggered();
  void on_unGroup_triggered();
  void on_func_test_triggered();
  void on_copy();
  void on_paste();
  void on_cut();
  void dataChanged();
  void positionChanged(int x, int y);

  void about();
  void itemClicked(QString s);

 protected:
  void closeEvent(QCloseEvent *event) Q_DECL_OVERRIDE;

 private:
  void createLayerWidgets();
  void createCommandWidget();
  void createActions();
  void createMenus();
  void createToolbars();
  void createPropertyEditor();
  void createToolBox();

  DrawView *activeMdiChild();
  QMdiSubWindow *findMdiChild(const QString &fileName);
  QTreeWidget *layerWidget;
  QListWidget *itemWidget;
  QMenu *windowMenu;
  QDockWidget *mapViewDock;
  QMdiArea *mdiArea;
  QSignalMapper *windowMapper;
  CommandWidget *mCommandWidget;
  // update ui
  QTimer m_timer;
  // toolbox
  QToolBox *toolBox;
  // edit toolbar;
  QToolBar *editToolBar;
  // align toolbar
  QToolBar *alignToolBar;
  //  QAction *rightAct;
  //  QAction *leftAct;
  //  QAction *vCenterAct;
  //  QAction *hCenterAct;
  //  QAction *upAct;
  //  QAction *downAct;
  //  QAction *horzAct;
  //  QAction *vertAct;
  //  QAction *heightAct;
  //  QAction *widthAct;
  //  QAction *allAct;
  //  QAction *bringToFrontAct;
  //  QAction *sendToBackAct;

  //  QAction *funcAct;
  // file
  QAction *newAct;
  QAction *openAct;
  QAction *saveAct;
  QAction *exitAct;

  QAction *groupAct;
  QAction *unGroupAct;

  // edit action
  QAction *deleteAct;
  QAction *undoAct;
  QAction *redoAct;
  QAction *copyAct;
  QAction *pasteAct;
  QAction *cutAct;
  QAction *zoomInAct;
  QAction *zoomOutAct;
  QAction *addlayer;
  QAction *removelayer;
  // drawing toolbar
  QToolBar *drawToolBar;
  QActionGroup *drawActionGroup;
  QAction *selectAct;
  QAction *dragAct;
  //    QAction  * lineAct;
  QAction *rectAct;
  //  QAction *roundRectAct;
  //  QAction *ellipseAct;
  QAction *polygonAct;
  QAction *polylineAct;
  QAction *textAct;
  //  QAction *bezierAct;
  //  QAction *rotateAct;

  QAction *closeAct;
  QAction *closeAllAct;
  QAction *tileAct;
  QAction *cascadeAct;
  QAction *nextAct;
  QAction *previousAct;
  QAction *separatorAct;
  QAction *aboutAct;
  QAction *aboutQtAct;

  // property editor
  QDockWidget *dockProperty;
  ObjectController *propertyEditor;
  QObject *theControlledObject;

  QListWidget *listView;

  QUndoStack *undoStack;
  QUndoView *undoView;
  // statusbar label
  QLabel *m_posInfo;
};

#endif  // MAINWINDOW_H
