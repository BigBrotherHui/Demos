#ifndef TITLEBAR_H
#define TITLEBAR_H

#include <QWidget>

namespace Ui {
class TitleBar;
}

class TitleBar : public QWidget
{
    Q_OBJECT

public:
    explicit TitleBar(QWidget *parent = nullptr);
    ~TitleBar();
    void setTitle(const QString &title);
    QString getTitle();
    void paintEvent(QPaintEvent *event) override;
    bool event(QEvent *event) override;
    void setTarget(QWidget *target);
    void hideCloseButton();
    void closeByOthers();
protected:
    void take_over();
    void slot_targetChanged(QWidget *w);
signals:
    void signal_propertyChanged(int value);
    void signal_targetChanged(QWidget *w);
    void signal_closeWindow();
public slots:
    void slot_propertyChanged(int value);
private:
    Ui::TitleBar *ui;
    QWidget *target{nullptr};
};

#endif // TITLEBAR_H
