#include "ufunbaseinit.h"
#include <QSettings>
#include <QDebug>
#include "ustatus.h"

uFunBaseInit *uFunBaseInit::self = nullptr;

uFunBaseInit::uFunBaseInit()
{
    ;
}

uFunBaseInit * uFunBaseInit::getInStance()
{
    if (self == nullptr)
    {
        if (self == nullptr)
        {
            self = new uFunBaseInit;
        }
    }
    return self;
}

void uFunBaseInit::f_InitBaseWidgetName()
{
    QString mIniFileName = uStatus::mMainPath + "/"  + uStatus::mMainConfigureFileName;
    QSettings *mMainConfigIni = new QSettings(mIniFileName, QSettings::IniFormat);
    if(!mMainConfigIni->contains("/base/mainwindow"))
        mMainConfigIni->setValue("/base/mainwindow","MainWindow");
    uStatus::mMainWindow = mMainConfigIni->value("/base/mainwindow").toString();
    if(!mMainConfigIni->contains("/base/MainToolbarWidget"))
        mMainConfigIni->setValue("/base/MainToolbarWidget","MainToolBar");
    uStatus::mMainToolbarWidget = mMainConfigIni->value("/base/MainToolbarWidget").toString();
    if(!mMainConfigIni->contains("/base/MainCenterWidget"))
        mMainConfigIni->setValue("/base/MainCenterWidget","CaseManage");
    uStatus::mMainCenterWidget = mMainConfigIni->value("/base/MainCenterWidget").toString();
    if(!mMainConfigIni->contains("/base/MainControlWidgetRight"))
        mMainConfigIni->setValue("/base/MainControlWidgetRight","CaseManageControl");
    uStatus::mMainControlWidgetRight = mMainConfigIni->value("/base/MainControlWidgetRight").toString();
    if(!mMainConfigIni->contains("/base/MainControlWidgetLeft"))
        mMainConfigIni->setValue("/base/MainControlWidgetLeft","NULL");
    uStatus::mMainControlWidgetLeft = mMainConfigIni->value("/base/MainControlWidgetLeft").toString();
    if(!mMainConfigIni->contains("/base/MainFootbarWidget"))
        mMainConfigIni->setValue("/base/MainFootbarWidget","NULL");
    uStatus::mMainFootbarWidget = mMainConfigIni->value("/base/MainFootbarWidget").toString();

    delete  mMainConfigIni;
}
