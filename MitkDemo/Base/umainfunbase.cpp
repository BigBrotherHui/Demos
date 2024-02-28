#include "umainfunbase.h"
#include <QPointer>
#include "ufunbase.h"
#include "ufunction.h"
#include "ustatus.h"
#include "ufunmap.h"
#include <QDebug>
uMainFunBase::uMainFunBase()
{

}

void uMainFunBase::f_SetStackedWidget_Control_right(QStackedWidget *pStackedWidget_control)
{
    mStackedWidget_control_right = pStackedWidget_control;
}
void uMainFunBase::f_SetStackedWidget_Control_left(QStackedWidget *pStackedWidget)
{
    mStackedWidget_control_left = pStackedWidget;
}
void uMainFunBase::f_SetStackedWidget_Control_toolbar(QStackedWidget *pStackedWidget_control)
{
    mStackedWidget_control_toolbar = pStackedWidget_control;
}
void uMainFunBase::f_SetStackedWidget_Control_footbar(QStackedWidget *pStackedWidget_control)
{
    mStackedWidget_control_footbar = pStackedWidget_control;
}
void uMainFunBase::f_SetStackedWidget_Center(QStackedWidget *pStackedWidget)
{
    mStackedWidget_Center = pStackedWidget;
}

QStackedWidget *uMainFunBase::f_GetStackedWidget_Control_right()
{
    return mStackedWidget_control_right;
}

QStackedWidget *uMainFunBase::f_GetStackedWidget_Control_left()
{
    return mStackedWidget_control_left;
}

QStackedWidget *uMainFunBase::f_GetStackedWidget_Control_toolbar()
{
    return mStackedWidget_control_toolbar;
}

QStackedWidget *uMainFunBase::f_GetStackedWidget_Control_footbar()
{
    return mStackedWidget_control_footbar;
}

QStackedWidget *uMainFunBase::f_GetStackedWidget_Center()
{
    return mStackedWidget_Center;
}




void uMainFunBase::f_Open_Center(QString pName)
{
    f_Open_Center(pName,"");
}
void uMainFunBase::f_Open_Center(QString pName,QString pParameter)
{
    if(pName=="NULL")
    {
        mStackedWidget_Center->setVisible(false);
        return;
    }
    mStackedWidget_Center->setVisible(true);
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance(pName+"*","");
    for (int i =0 ; i < mStackedWidget_Center->count() ; i++)
    {
        if (mStackedWidget_Center->widget(i) == mCur)
        {
            mStackedWidget_Center->setCurrentWidget(mCur);
            return;
        }
    }

//    mCur->f_Init();
    if (!pParameter.isEmpty())
    {
        mCur->f_SetPparameter(pParameter);
    }

    mStackedWidget_Center->addWidget( mCur );
    mStackedWidget_Center->setCurrentWidget(mCur);
}



void uMainFunBase::f_OpenControl_right(QString pName)
{
    f_OpenControl_right(pName,"");
}
void uMainFunBase::f_OpenControl_right(QString pName,QString pParameter)
{
    if(pName=="NULL")
    {
        mStackedWidget_control_right->setVisible(false);
        return;
    }
    mStackedWidget_control_right->setVisible(true);
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance(pName+"*",pParameter);

    for (int i =0 ; i < mStackedWidget_control_right->count() ; i++)
    {
        if (mStackedWidget_control_right->widget(i) == mCur)
        {
            mCur->f_Refresh();
            mStackedWidget_control_right->setCurrentWidget(mCur);
            return;
        }
    }
    mCur->f_Init();
    if (!pParameter.isEmpty())
    {
        mCur->f_SetPparameter(pParameter);
    }
    mStackedWidget_control_right->addWidget( mCur );
    mStackedWidget_control_right->setCurrentWidget(mCur);
}

void uMainFunBase::f_OpenControl_left(QString pName)
{
    f_OpenControl_left(pName,"");
}
void uMainFunBase::f_OpenControl_left(QString pName,QString pParameter)
{
    if(pName=="NULL")
    {
        mStackedWidget_control_left->setVisible(false);
        return;
    }
    mStackedWidget_control_left->setVisible(true);
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance(pName+"*",pParameter);

    for (int i =0 ; i < mStackedWidget_control_left->count() ; i++)
    {
        if (mStackedWidget_control_left->widget(i) == mCur)
        {

            mStackedWidget_control_left->setCurrentWidget(mCur);
            return;
        }
    }

//    mCur->f_Init();
    if (!pParameter.isEmpty())
    {
        mCur->f_SetPparameter(pParameter);
    }

    mStackedWidget_control_left->addWidget( mCur );
    mStackedWidget_control_left->setCurrentWidget(mCur);
}

void uMainFunBase::f_OpenControl_toobar(QString pName)
{
    f_OpenControl_toobar(pName,"");
}
void uMainFunBase::f_OpenControl_toobar(QString pName,QString pParameter)
{
    if(pName=="NULL")
    {
        mStackedWidget_control_toolbar->setVisible(false);
        return;
    }
    mStackedWidget_control_toolbar->setVisible(true);
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance(pName+"*",pParameter);
    for (int i =0 ; i < mStackedWidget_control_toolbar->count() ; i++)
    {
        if (mStackedWidget_control_toolbar->widget(i) == mCur)
        {

            mStackedWidget_control_toolbar->setCurrentWidget(mCur);
            return;
        }
    }

//    mCur->f_Init();
    if (!pParameter.isEmpty())
    {
        mCur->f_SetPparameter(pParameter);
    }

    mStackedWidget_control_toolbar->addWidget( mCur );
    mStackedWidget_control_toolbar->setCurrentWidget(mCur);
}

void uMainFunBase::f_OpenControl_footbar(QString pName)
{
    f_OpenControl_footbar(pName,"");
}
void uMainFunBase::f_OpenControl_footbar(QString pName,QString pParameter)
{
    if(pName=="NULL")
    {
        mStackedWidget_control_footbar->setVisible(false);
        return;
    }
    mStackedWidget_control_footbar->setVisible(true);
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance(pName+"*",pParameter);

    for (int i =0 ; i < mStackedWidget_control_footbar->count() ; i++)
    {
        if (mStackedWidget_control_footbar->widget(i) == mCur)
        {

            mStackedWidget_control_footbar->setCurrentWidget(mCur);
            return;
        }
    }

//    mCur->f_Init();
    if (!pParameter.isEmpty())
    {
        mCur->f_SetPparameter(pParameter);
    }

    mStackedWidget_control_footbar->addWidget( mCur );
    mStackedWidget_control_footbar->setCurrentWidget(mCur);
}


void uMainFunBase::f_CloseControl_right(QString pName)
{
    f_CloseControl_right(pName,"");
}
void uMainFunBase::f_CloseControl_right(QString pName,QString pParameter)
{
    if(pName=="NULL")
    {
        return;
    }

    QMap<QString,QPointer<uFunBase>>::Iterator it= uFunMap::mFunBaseMap->find(pName+pParameter);

    if(it!=uFunMap::mFunBaseMap->end())
    {
        qDebug() << "uFunction::f_GetObjectInstance Exist";

        QPointer<uFunBase> an1 = (QPointer<uFunBase>)(it.value());
        if (an1)
        {
            mStackedWidget_Center->removeWidget( an1 );
        }

    }
    uFunction::getInStance()->f_RemveObjectInstance(pName+"*",pParameter);

}


void uMainFunBase::f_CloseControl_Center(QString pName)
{
    f_CloseControl_Center(pName,"");
}
void uMainFunBase::f_CloseControl_Center(QString pName,QString pParameter)
{
    if(pName=="NULL")
    {
        return;
    }

    QMap<QString,QPointer<uFunBase>>::Iterator it= uFunMap::mFunBaseMap->find(pName+pParameter);

    if(it!=uFunMap::mFunBaseMap->end())
    {
        qDebug() << "uFunction::f_GetObjectInstance Exist";

        QPointer<uFunBase> an1 = (QPointer<uFunBase>)(it.value());
        if (an1)
        {
            mStackedWidget_Center->removeWidget( an1 );
        }

    }

    uFunction::getInStance()->f_RemveObjectInstance(pName+"*",pParameter);
}

void uMainFunBase::f_InitControl(QString pName)
{
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance(pName+"*","");
}
