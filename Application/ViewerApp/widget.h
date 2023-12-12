#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include <QMenu>

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
    void switchWidget(QString pluginName);
private slots:
    void on_pushButton_clicked();
    void on_listWidget_customContextMenuRequested(const QPoint& pos);
    void on_pushButton_3_clicked();
    void on_pushButton_broadcast_clicked();
protected slots:
    void slot_unloadAction_triggered();
    void slot_listWidget_clicked(const QModelIndex& index);
private:
    Ui::Widget *ui;
    QMenu* m_popMenu;
};
#endif // WIDGET_H
