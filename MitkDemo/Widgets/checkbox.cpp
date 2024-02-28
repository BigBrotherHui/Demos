#include "checkbox.h"
#include <QDebug>
#include <QPaintEvent>
#include <QPainter>
#include <QPropertyAnimation>
#include <QButtonGroup>
CheckBox::CheckBox(QWidget *w) : QCheckBox(w), pose(startvalue) {
  animation = new QPropertyAnimation();
  animation->setPropertyName("m_pos");
  animation->setTargetObject(this);
  animation->setEasingCurve(QEasingCurve::Linear);
  animation->setDuration(100);
  animation->setLoopCount(1);
  take_over();
}

CheckBox::~CheckBox() {
  if (animation) {
    delete animation;
    animation = nullptr;
  }
  setFixedSize(48,24);
}

void CheckBox::paintEvent(QPaintEvent *) {
  bool ischecked = isChecked();
  QPainter painter(this);
  painter.setRenderHints(QPainter::Antialiasing);
  painter.setViewport(rect());
  painter.setWindow(0, 0, 48, 24);
  QColor color,colorBrush;
  if (!ischecked)
  {
    colorBrush = handlecolorUnchecked;
    color=bkcolorUnchecked;
  }
  else
  {
    colorBrush = handlecolorChecked;
    color=bkcolorChecked;
  }
  painter.setBrush(QBrush(color));
  painter.setPen(color);
  painter.drawRoundedRect(0, 0, 48, 24, 12, 12);
  painter.setBrush(colorBrush);
  if(animation && animation->state()!=QPropertyAnimation::Running)
  {
      if(ischecked)
      {
          pose=endvalue;
      }
      else
      {
          pose=startvalue;
      }
  }
  painter.drawEllipse(pose, 2, 20, 20);
}

int CheckBox::getPos() { return pose; }

void CheckBox::setPos(int pos)
{
  this->pose = pos;
  update();
}

void CheckBox::clicked() {
    if (animation->state() == QPropertyAnimation::Running || (isExclusive && isChecked()))
    {
            return;
    }
    if (isChecked()) {
        animation->setStartValue(endvalue);
        animation->setEndValue(startvalue);
    } else {
        animation->setStartValue(startvalue);
        animation->setEndValue(endvalue);
    }
    animation->start();
    setChecked(!isChecked());
    emit QCheckBox::clicked(isChecked());
}

void CheckBox::SetCheckedBkColor(QColor color)
{
    bkcolorChecked=color;
}

void CheckBox::SetCheckedHandleColor(QColor color)
{
    handlecolorChecked=color;
}

void CheckBox::SetUncheckedBkColor(QColor color)
{
    bkcolorUnchecked=color;
}

void CheckBox::SetUncheckedHandleColor(QColor color)
{
    handlecolorUnchecked=color;
}

void CheckBox::SetExclusive(bool b)
{
    isExclusive=b;
}

bool CheckBox::event(QEvent *e)
{
    if (e->type() == QEvent::MouseButtonPress || e->type() == QEvent::MouseButtonDblClick)
    {
        QMouseEvent *me=static_cast<QMouseEvent *>(e);
        if(me->button() != Qt::LeftButton )
            return false;
        clicked();
        return true;
    }
    return QCheckBox::event(e);
}

void CheckBox::take_over() {}
