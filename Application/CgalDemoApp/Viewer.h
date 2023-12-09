#ifndef VIEWER_H
#define VIEWER_H
#include <QMap>
#include <CGAL/Qt/qglviewer.h>
#include "Scene.h"

// forward declarations
class QWidget;
class Viewer : public CGAL::QGLViewer{

  Q_OBJECT

public:
  Viewer(QWidget * parent);

  // overload several CGAL::QGLViewer virtual functions
  void draw();
  void initializeGL();
  void setScene(Scene* pScene);
protected:
  virtual void mousePressEvent(QMouseEvent* e);
  virtual void mouseReleaseEvent(QMouseEvent* e);

private:
  bool m_custom_mouse;
  Scene* m_pScene;
}; // end class Viewer

#endif // VIEWER_H
