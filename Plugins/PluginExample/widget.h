#ifndef WIDGET_H
#define WIDGET_H

#include "PluginInterface.h"

class PluginManager;
QT_BEGIN_NAMESPACE
namespace Ui { class Widget; }
QT_END_NAMESPACE

class Widget : public WidgetBase
{
    Q_OBJECT

public:
    Widget(QWidget *parent = nullptr);
    ~Widget();
protected:
private slots:

private:
    Ui::Widget *ui;
};
#endif // WIDGET_H
