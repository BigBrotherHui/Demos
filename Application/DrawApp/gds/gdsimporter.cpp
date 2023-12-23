#include "gdsimporter.h"
#include <fstream>
#include <sstream>
#include "graphicsrectitem.h"
#include "drawscene.h"
#include "gdsCpp.hpp"
#include "graphicsbezier.h"
#include <QDebug>
gdsimporter::gdsimporter() {}

gdsimporter *gdsimporter::instance() {
  static gdsimporter s;
  return &s;
}
QRect createQRect(const std::vector<int> &points) {
  if (points.size() >= 8) {
    int x = points[0];
    int y = points[1];
    int width = points[2] - points[0];
    int height = points[5] - points[1];
    return QRect(x, y, width, height);
  }

  return QRect();
}
void gdsimporter::parseFile(QString path, std::vector<GraphicsItem *> &items) {
  std::ifstream file(path.toStdString());
  std::string line;
  std::vector<int> points;
  while (std::getline(file, line)) {
    if (line.substr(0, 5) == "[XY]:") {
      size_t start = line.find('{') + 1;
      size_t end = line.find('}');
      std::string data = line.substr(start, end - start);
      std::istringstream iss(data);
      int value;
      points.clear();
      while (iss >> value) {
        points.push_back(value);
        char comma;
        iss >> comma;  // 读取逗号
      }
      QRect rect = createQRect(points);
      GraphicsRectItem *r = new GraphicsRectItem(rect);
      items.push_back(r);
    }
  }
  file.close();
}

void gdsimporter::writeFile(QString path, const std::vector<GraphicsItem *> &items)
{
    gdscpp writer;
    gdsSTR str;
    for(int i=0;i<items.size();++i)
    {
        GraphicsItem *item=items[i];
        switch( item->itemtype())
        {
        case rectangle:
        {
            GraphicsRectItem *rectitem=static_cast<GraphicsRectItem *>(item);
            str.name="box";
            str.BOUNDARY.push_back(draw2ptBox(rectitem->zValue(),rectitem->pos().x(),rectitem->pos().y(),rectitem->pos().x()+rectitem->width(),
                                              rectitem->pos().y()+rectitem->height()));
            break;
        }
        case polygon: {
            GraphicsPolygonItem *polygonitem=static_cast<GraphicsPolygonItem *>(item);
            str.name="polygon";
            std::vector<int> x,y;
            std::vector<QPointF> pts=polygonitem->points().toStdVector();
            for(int k=0;k<pts.size();++k)
            {
                x[i]=pts[i].x();
                y[i]=pts[i].y();
            }
            str.BOUNDARY.push_back(drawBoundary(polygonitem->zValue(),x,y));
            break;
        }
        case polyline: {
            GraphicsBezier *pathitem=static_cast<GraphicsBezier *>(item);
            str.name="polyline";
            std::vector<int> x,y;
            std::vector<QPointF> pts=pathitem->points().toStdVector();
            for(int k=0;k<pts.size();++k)
            {
                x[i]=pts[i].x();
                y[i]=pts[i].y();
            }
            str.BOUNDARY.push_back(drawBoundary(pathitem->zValue(),x,y));
            break;
        }
        case text:
        {

            break;
        }
        case DrawShape::instance:
        {

            break;
        }
        default:
        break;
        }
    }
    writer.setSTR(str);
    writer.write(path.toStdString());
}
