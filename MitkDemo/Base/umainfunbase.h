#ifndef UMAINFUNBASE_H
#define UMAINFUNBASE_H

#include <QString>
#include <QStackedWidget>
#ifdef QT_BOOTSTRAPPED
#  define QT_NO_EXCEPTIONS
#endif
#if !defined(QT_NO_EXCEPTIONS) && defined(Q_CC_GNU) && !defined (__EXCEPTIONS) && !defined(Q_MOC_RUN)
#  define QT_NO_EXCEPTIONS
#endif

#ifdef QT_NO_EXCEPTIONS
#  define QT_TRY if (true)
#  define QT_CATCH(A) else
#  define QT_THROW(A) qt_noop()
#  define QT_RETHROW qt_noop()
#else
#  define QT_TRY try
#  define QT_CATCH(A) catch (A)
#  define QT_THROW(A) throw A
#  define QT_RETHROW throw
#endif

class uMainFunBase
{
public:
    uMainFunBase();

    void f_SetStackedWidget_Control_right(QStackedWidget *pStackedWidget_control);
    void f_SetStackedWidget_Control_left(QStackedWidget *pStackedWidget_control);
    void f_SetStackedWidget_Control_toolbar(QStackedWidget *pStackedWidget_control);
    void f_SetStackedWidget_Control_footbar(QStackedWidget *pStackedWidget_control);
    void f_SetStackedWidget_Center(QStackedWidget *pStackedWidget);

    QStackedWidget *f_GetStackedWidget_Control_right();
    QStackedWidget *f_GetStackedWidget_Control_left();
    QStackedWidget *f_GetStackedWidget_Control_toolbar();
    QStackedWidget *f_GetStackedWidget_Control_footbar();
    QStackedWidget *f_GetStackedWidget_Center();

    virtual void f_Open_Center(QString pName);
    virtual void f_OpenControl_left(QString pName);
    virtual void f_OpenControl_right(QString pName);
    virtual void f_OpenControl_toobar(QString pName);
    virtual void f_OpenControl_footbar(QString pName);
    virtual void f_Open_Center(QString pName,QString pParameter);
    virtual void f_OpenControl_left(QString pName,QString pParameter);
    virtual void f_OpenControl_right(QString pName,QString pParameter);
    virtual void f_OpenControl_toobar(QString pName,QString pParameter);
    virtual void f_OpenControl_footbar(QString pName,QString pParameter);

    virtual void f_CloseControl_right(QString pName);
    virtual void f_CloseControl_right(QString pName,QString pParameter);
    virtual void f_CloseControl_Center(QString pName);
    virtual void f_CloseControl_Center(QString pName,QString pParameter);

    virtual void f_InitControl(QString pName);

protected:

    QStackedWidget *mStackedWidget_control_left;
    QStackedWidget *mStackedWidget_control_right;
    QStackedWidget *mStackedWidget_control_toolbar;
    QStackedWidget *mStackedWidget_control_footbar;
    QStackedWidget *mStackedWidget_Center;
};

#endif // UMAINFUNBASE_H
