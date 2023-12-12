#pragma once

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
protected slots:
    void slot_simplify_clicked();
    void slot_cut_clicked();
    void slot_fillhole_clicked();
    void slot_load_clicked();
    void slot_offset_clicked();
private slots:

private:
    Ui::Widget *ui;
    class Impl;
    std::unique_ptr<Impl> m_impl;
};
