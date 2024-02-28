#ifndef SYSTEMTOOLBAR_H
#define SYSTEMTOOLBAR_H

#include <QWidget>

namespace Ui {
class SystemToolBar;
}

class SystemToolBar : public QWidget
{
    Q_OBJECT

public:
    explicit SystemToolBar(QWidget *parent = nullptr);
    ~SystemToolBar();

private slots:
    void on_pushButton_settings_clicked();

    void on_pushButton_screenshot_clicked();

    void on_pushButton_exit_clicked();

private:
    Ui::SystemToolBar *ui;
};

#endif // SYSTEMTOOLBAR_H
