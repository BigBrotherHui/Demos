#include <QDebug>
#include <QPointer>
#include <QMainWindow>
#include "ufunction.h"
#include "ustatus.h"
#include "ufunbase.h"
#include "ufunmap.h"

uFunction::uFunction()
{

}


uFunction *uFunction::self = nullptr;
uFunction * uFunction::getInStance()
{
    if (self == nullptr)
    {
        if (self == nullptr) {
            self = new uFunction;
        }
    }
    return self;
}


void uFunction::f_Test()
{
    qDebug() << "uFunction::f_Test ... ...";
}

uMainFunBase *uFunction::f_GetMain()
{
    return mMain;
}

void uFunction::f_InsertMain(uMainFunBase *pMain)
{
    mMain = pMain;
}

QPointer<uFunBase> uFunction::f_GetObjectInstance(QString pObjectName, QString pParameter)
{
    qDebug() << "uFunction::f_GetObjectInstance ... ... : " << pObjectName << "," << pParameter;
    QMap<QString,QPointer<uFunBase>>::Iterator it= uFunMap::mFunBaseMap->find(pObjectName+pParameter);

    if(it!=uFunMap::mFunBaseMap->end())
    {
        qDebug() << "uFunction::f_GetObjectInstance Exist";

        QPointer<uFunBase> an1 = (QPointer<uFunBase>)(it.value());
        return an1;
    }
    else
    {
        qDebug() << "uFunction::f_GetObjectInstance not Exist";
        int type = QMetaType::type( (pObjectName ).toStdString().c_str());
        qDebug() << type;
        const QMetaObject *metaObj = QMetaType::metaObjectForType(type);
        QObject *obj = metaObj->newInstance(Q_ARG(QWidget*, nullptr));
        qDebug() <<"uFunction::f_GetObjectInstance "<< obj;

        QPointer<uFunBase> o_an = QPointer<uFunBase>(qobject_cast<uFunBase*>(obj));
        qDebug() <<"uFunction::f_GetObjectInstance Insert :"<< pObjectName+pParameter;
        uFunMap::mFunBaseMap->insert(pObjectName+pParameter , o_an);

//        uFunBase *an = qobject_cast<uFunBase*>(obj);
        o_an->f_Init();
        return o_an;
    }
}


bool uFunction::f_RemveObjectInstance(QString pObjectName,QString pParameter)
{
    qDebug() << "uFunction::f_RemveObjectInstance ... ... : " << pObjectName << "," << pParameter;

    QMap<QString,QPointer<uFunBase>>::Iterator it= uFunMap::mFunBaseMap->find(pObjectName+pParameter);
    if(it==uFunMap::mFunBaseMap->end())
    {
        qDebug() << "uFunction::f_RemveObjectInstance Not Exist";
        return true;
    }
    qDebug() << "uFunction::f_RemveObjectInstance Exist , Key : " << it.key();
    QPointer<uFunBase> an1 = (QPointer<uFunBase>)(it.value());
    qDebug() << "uFunction::f_RemveObjectInstance Exist , Key : " << it.key()<<"f_BeforeClose";
    an1->f_BeforeClose();
    an1->f_Close();
    an1->deleteLater();
    uFunMap::mFunBaseMap->remove(it.key());
}
