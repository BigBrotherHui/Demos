#ifndef RULER_H
#define RULER_H

#include <QWidget>

namespace Ui {
class Ruler;
}

class Ruler : public QWidget {
  Q_OBJECT

 public:
  explicit Ruler(QWidget *parent = nullptr);
  ~Ruler();
  void updateValue(int v);
 signals:
  void signal_zoomFactorChanged(int);

 private:
  Ui::Ruler *ui;
};

#endif  // RULER_H
