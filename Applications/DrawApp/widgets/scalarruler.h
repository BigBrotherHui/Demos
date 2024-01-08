#ifndef SCALARRULER_H
#define SCALARRULER_H

#include <QtWidgets>
//#include <QWidget>
//#include <QPainter>
#include <QAbstractScrollArea>
#include <QResizeEvent>
#include <QScrollBar>

class ScalarRuler : public QWidget {
  Q_OBJECT
  QScrollBar* scroll;

 public:
  //    Ruler(QAbstractScrollArea* parent=nullptr): QWidget(parent),offset(0)
  ScalarRuler(QAbstractScrollArea* parent = nullptr)
      : QWidget(parent), scroll(new QScrollBar(this)) {
    setFixedSize(40, height());
    move(0, 40);
    connect(scroll, &QScrollBar::valueChanged, this, &ScalarRuler::setOffset);
    // connect(parent->verticalScrollBar(), &QScrollBar::valueChanged, this,
    // &Ruler::setOffset); connect(this, SIGNAL(valueChanged(int)), this,
    // SLOT(setOffset(int)));
  }
  virtual void paintEvent(QPaintEvent* event) {
    QPainter painter(this);
    painter.translate(0, -offset);
    int const heightMM = height() * toMM();
    painter.setFont(font());
    QFontMetrics fm(font());
    for (int position = 0; position < heightMM; ++position) {
      int const positionInPix = int(position / toMM());
      if (position % 10 == 0) {
        if (position != 0) {
          QString const txt = QString::number(position);
          QRect txtRect = fm.boundingRect(txt).translated(0, positionInPix);
          txtRect.translate(0, txtRect.height() / 2);
          painter.drawText(txtRect, txt);
        }
        painter.drawLine(width() - 15, positionInPix, width(), positionInPix);
      } else {
        painter.drawLine(width() - 10, positionInPix, width(), positionInPix);
      }
    }
  }

  virtual void resizeEvent(QResizeEvent* event) {
    int const maximumMM = event->size().height() * toMM();
    QFontMetrics fm(font());
    int w = fm.width(QString::number(maximumMM)) + 20;
    if (w != event->size().width()) {
      QSize const newSize(w, event->size().height());
      sizeChanged(newSize);
      return setFixedSize(newSize);
    }
    return QWidget::resizeEvent(event);
  }

  void setOffset(int value) {
    offset = value;
    update();
  }
 signals:
  void sizeChanged(QSize const&);

 private:
  int offset;

  static qreal toMM() { return 25.4 / qApp->desktop()->logicalDpiY(); }
};
#endif  // RULER_H
