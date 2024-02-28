#pragma execution_character_set("utf-8")
#include "login.h"
#include "ui_login.h"

// qt
#include <QCryptographicHash>
#include <QDebug>
#include <QKeyEvent>
#include <QRegExp>
#include <QRegExpValidator>
#include <QSettings>
#include <QSharedPointer>
#include <QSqlDatabase>
#include <QSqlError>
#include <QSqlQuery>
#include <QMessageBox>
#include <QTimer>
#include <QDateTime>
// vtk
#include <itkOutputWindow.h>
#include <vtkOutputWindow.h>

#include "runtimeutil.h"
#include "systemutil.h"
#include "SqlDemo.h"
#include "global.h"
#include "ustatus.h"
#include "ufunction.h"
#include "mainframe.h"
#include "cryptoutil.h"
Login::Login(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Login)
{
    ui->setupUi(this);
    if(!mSql)
    {
        mSql=new SqlDemo;
    }
    qDebug() << "init login page...";

    SystemUtil::InstallEventFitler(this);
    if (RuntimeUtil::GetRememberMe())
    {
          ui->checkBox_rememberMe->setCheckState(Qt::Checked);
          ui->lineEdit_loginName->setText(RuntimeUtil::GetLoginName());
          ui->lineEdit_password->setText("password");
          mPassword = RuntimeUtil::GetPassword();
    }
    else if (RuntimeUtil::GetDebug())
    {
        ui->lineEdit_loginName->setText(LOGIN_ID_ADMIN);
        ui->lineEdit_password->setText(DEFAULT_PASSWORD);
    }
    //设置账号、密码框图标
    {
        QHBoxLayout *layout = new QHBoxLayout(this);
        QLabel *label = new QLabel(this);
        layout->addWidget(label);
        layout->setContentsMargins(28, 0, 0, 0);
        ui->lineEdit_loginName->setLayout(layout);
        ui->lineEdit_loginName->setTextMargins(58, 0, 0, 0);
        ui->lineEdit_loginName->setLayout(layout);
    }
    {
        QHBoxLayout *layout = new QHBoxLayout(this);
        QLabel *label = new QLabel(this);
        layout->addWidget(label);
        layout->setContentsMargins(28, 0, 0, 0);
        ui->lineEdit_password->setLayout(layout);
        ui->lineEdit_password->setTextMargins(58, 0, 0, 0);
        ui->lineEdit_password->setLayout(layout);
    }
    {
        //初始化closeButton
//        ui->closeButton->setIconDirection(IconButton::RIGHT);
    }
    ui->lineEdit_loginName->setFocus();
    ui->closeButton->setMinimumSize(96,34);
}

Login::~Login()
{
    qDebug() << "Login::~Login";
    delete ui;
    if(mSql)
    {
        delete mSql;
        mSql=nullptr;
    }
    qDebug() << "Login::~Login finished";
}

bool Login::eventFilter(QObject *target, QEvent *event)
{
    if (event->type() == QEvent::KeyPress)
    {
        QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
        if (keyEvent->key() == Qt::Key_Return || keyEvent->key() == Qt::Key_Enter)
        {
            // 回车处理：切换焦点，触发登录
            if (target == ui->lineEdit_loginName)
            {
                ui->lineEdit_password->setFocus();
            }
            else if (target == ui->lineEdit_password||target==ui->pushButton_login)
            {
                on_pushButton_login_clicked();
            }
            return true;
        }
        else if (keyEvent->key() == Qt::Key_Escape)
        {
            // ESC键处理：退出系统
            on_closeButton_clicked();
            return true;
        }
      }
    else if (event->type() == QEvent::FocusIn)
    {
        if (target->inherits("QLineEdit"))
        {
            QLineEdit *edit = static_cast<QLineEdit *>(target);
            edit->selectAll();
        }
    }
    return QObject::eventFilter(target, event);
}

void Login::on_pushButton_login_clicked()
{
    qDebug() << "Login button clicked...";
    QString loginName = ui->lineEdit_loginName->text();
    QString password = ui->lineEdit_password->text();

    if (loginName.isEmpty())
    {
        QMessageBox::warning(this, tr("错误"), tr("请输入用户名！"));
        ui->lineEdit_loginName->setFocus();
        return;
    }

    if ( loginName.length()>20 )
    {
        QMessageBox::warning(this, tr("错误"), tr("用户名长度必须在20位以内！"));
        ui->lineEdit_loginName->setFocus();
        return;
    }
    if (password.isEmpty())
    {
        QMessageBox::warning(this, tr("错误"), tr("请输入密码！"));
        ui->lineEdit_password->setFocus();
        return;
    }

    QRegExp regx(QString("[%1]+").arg(gSpecial));
    if (regx.indexIn(loginName) != -1 || regx.indexIn(password) != -1)
    {
        return;
    }
    QString passwordMd5 = mPassword.isEmpty()?CryptoUtil::md5(password):mPassword;
    bool isExisted{false};
    bool ret=mSql->IsUserExisted(loginName,isExisted);
    if(!ret)
    {
        QMessageBox::warning(this,"提示","查询用户是否存在失败！");
        return;
    }
    if (!isExisted)
    {
        QMessageBox::warning(this, tr("错误"), tr("无此账号！"));
        return;
    }
    int enable;
    ret=mSql->IsUserEnable(loginName,enable);
    if(!ret)
    {
        QMessageBox::warning(nullptr,"提示","查询是否可用失败！");
        return;
    }
    else
    {
        if(enable==Enable_Disable)
        {
            QMessageBox::warning(nullptr,"提示","该账号已被停用！");
            return;
        }
    }
    user_info info;
    ret=mSql->SelectUser(loginName,passwordMd5,info);
    if(!ret)
    {
        QMessageBox::warning(this, tr("错误"), tr("密码错误！"));
        return;
    }
    gCurrentUser=info;
    bool rememberMe = ui->checkBox_rememberMe->isChecked();
    RuntimeUtil::SetRememberMe(rememberMe);
    if (rememberMe)
    {
        RuntimeUtil::SetLastLogin(info.id, info.password);
    }
    emit signal_loginSucceed();
    return;
}

void Login::on_closeButton_clicked()
{
    qDebug() << "exit button down...";
    if (QMessageBox::Yes == QMessageBox::question(this, tr("提示"), tr("确认要退出系统？"), QMessageBox::Yes | QMessageBox::No))
    {
        uStatus::mSql->DisConnectDB();
        qApp->quit();
    }
}

void Login::on_lineEdit_password_textChanged(const QString &pText)
{
    mPassword = "";
    QLineEdit *lineedit=static_cast<QLineEdit *>(sender());
    QFont font=lineedit->font();
    if(pText.length()==0)
    {
        lineedit->setStyleSheet("font-size:24px;");
        font.setLetterSpacing(QFont::AbsoluteSpacing,0);
    }
    else
    {
        lineedit->setStyleSheet("font-size:15px;");
        font.setLetterSpacing(QFont::AbsoluteSpacing,10);
    }
    lineedit->setFont(font);
}
