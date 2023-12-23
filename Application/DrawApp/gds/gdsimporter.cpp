#include "gdsimporter.h"
#include <fstream>
#include <sstream>
#include "graphicsrectitem.h"
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
void gdsimporter::parseFile(QString path, std::vector<QGraphicsItem *> &items) {
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
