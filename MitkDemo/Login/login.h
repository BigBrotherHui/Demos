#ifndef LOGIN_H
#define LOGIN_H

#include <QDialog>
#include <memory>
#include <ufunbase.h>
using namespace std;

QT_BEGIN_NAMESPACE
namespace Ui {
class Login;
}
QT_END_NAMESPACE
class SqlDemo;
class MainWindow;

class Login : public QDialog
{
    Q_OBJECT

public:
    Login(QWidget *parent = nullptr);
    ~Login();
signals:
    void signal_loginSucceed();
private slots:
    /**
     * @brief 事件过滤器，用于屏蔽PlainTextEdit中的换行，以及增加回车切换焦点功能
     * @param target 事件目标对象
     * @param event 事件
     * @return true:过滤事件处理
     */
    bool eventFilter(QObject *target, QEvent *event);

    void on_pushButton_login_clicked();

    void on_closeButton_clicked();

    void on_lineEdit_password_textChanged(const QString &pText);

private:
    Ui::Login *ui;
    QString mPassword;
    SqlDemo *mSql{nullptr};
};
#endif  // LOGIN_H
