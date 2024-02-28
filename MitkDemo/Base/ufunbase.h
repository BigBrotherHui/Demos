#ifndef UFUNBASE_H
#define UFUNBASE_H
#include <QWidget>

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

class uFunBase: public QWidget
{
    Q_OBJECT
public:
    explicit uFunBase(QWidget *parent = 0);
    QString mID = "";
    QString mInstensName = "";
    QString mType="C";
    //
    QString f_GetType();
    void f_SetType(QString pType);


    virtual int f_BeforeClose();

    //
    virtual QString f_GetID();
    virtual void f_Init();
    virtual void f_Update();
    virtual void f_Destroy();
    virtual void f_Close();
    virtual void f_Refresh();
    virtual void f_Switch();
    virtual void f_SetPparameter(QString pParameter);
    void f_SetInstensName(QString pInstensName);
    QString f_GetInstensName();
};

#endif // UFUNBASE_H
