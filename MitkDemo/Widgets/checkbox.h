#ifndef CHECKBOX_H
#define CHECKBOX_H

#include <QCheckBox>
class QPropertyAnimation;
class CheckBox : public QCheckBox {
  Q_OBJECT
  Q_PROPERTY(int m_pos READ getPos WRITE setPos)
 public:
  CheckBox(QWidget *w);
  ~CheckBox();
  void paintEvent(QPaintEvent *) override;
  int getPos();
  void setPos(int pos);
  void clicked();
  void SetCheckedBkColor(QColor color);
  void SetCheckedHandleColor(QColor color);
  void SetUncheckedBkColor(QColor color);
  void SetUncheckedHandleColor(QColor color);
  void SetExclusive(bool b);
 protected:
  bool event(QEvent *e) override;
  void take_over();

 private:
  QPropertyAnimation *animation{nullptr};
  int pose{1};
  const static int startvalue{1};
  const static int endvalue{27};
  QColor handlecolorChecked{Qt::white},bkcolorChecked{"#1269D3"};
  QColor handlecolorUnchecked{QColor("#80838B")},bkcolorUnchecked{QColor("#2C313E")};
  QButtonGroup *g{nullptr};
  bool isExclusive{false};
};

#endif  // CHECKBOX_H
