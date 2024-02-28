#include "ufunbase.h"
#include <QUuid>
#include <QDebug>


uFunBase::uFunBase(QWidget *parent): QWidget(parent)
{
    QUuid id = QUuid::createUuid();
    uFunBase::mID = id.toString();
    qDebug() << "uFunBase UID :" << id.toString();
}
QString uFunBase::f_GetID()
{
    return mID;
}
void uFunBase::f_Init()
{

}
void uFunBase::f_Update()
{

}

QString uFunBase::f_GetType()
{
    return mType;
}

void uFunBase::f_SetType(QString pType)
{
    mType = pType;
}


void uFunBase::f_Switch()
{
    f_Update();
}

void uFunBase::f_Destroy()
{
    qDebug() << "uFunBase : f_Destroy ";
}

int uFunBase::f_BeforeClose()
{
    qDebug() << "uFunBase : f_BeforeClose ";
    return 0;
}

void uFunBase::f_Close()
{
    qDebug() << "uFunBase : f_Close ";

}

void uFunBase::f_Refresh()
{
    qDebug() << "uFunBase : f_Refresh ";
}


void uFunBase::f_SetPparameter(QString pParameter)
{
    qDebug() << "uFunBase : f_SePparameter ";
}

void uFunBase::f_SetInstensName(QString pInstensName)
{
    mInstensName = pInstensName;
}

QString uFunBase::f_GetInstensName()
{
    return mInstensName;
}
