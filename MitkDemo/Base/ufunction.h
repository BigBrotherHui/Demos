#ifndef UFUNCTION_H
#define UFUNCTION_H
#include <QObject>
#include "ufunbase.h"

class uMainFunBase;

class uFunction
{
public:
    static uFunction *getInStance();
    uFunction();
    uMainFunBase *mMain;
    void f_Test();
    QPointer<uFunBase> f_GetObjectInstance(QString pObjectName,QString pParameter);
    bool f_RemveObjectInstance(QString pObjectName,QString pParameter);
    void f_InsertMain(uMainFunBase *pMain);
    uMainFunBase *f_GetMain();


private:
    static uFunction *self;
};

#endif // UFUNCTION_H
