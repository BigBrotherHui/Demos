#pragma execution_character_set("utf-8")
#ifndef UFUNMAP_H
#define UFUNMAP_H
#include <QObject>
#include <QString>
#include <QPointer>

class uFunBase;
class uFunMap
{
public:
    uFunMap(){};

    /////////////////////////所有初始化的类///////////////////////////////////////
    static QMap<QString, QPointer<uFunBase> > *mFunBaseMap;
};

#endif // UFUNMAP_H
