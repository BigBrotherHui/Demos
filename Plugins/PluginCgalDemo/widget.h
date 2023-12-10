#pragma once

#include <QWidget>

class PluginManager;
QT_BEGIN_NAMESPACE
namespace Ui { class Widget; }
QT_END_NAMESPACE

class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent = nullptr);
    ~Widget();
protected slots:
    void slot_btn_clicked();
    void slot_cut_clicked();
    void slot_fillhole_clicked();
private slots:

private:
    Ui::Widget *ui;
    class Impl;
    std::unique_ptr<Impl> m_impl;
};
