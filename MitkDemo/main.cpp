#include <QApplication>
#include <QDesktopWidget>
#include <QDebug>
#include <QPointer>
#include <QSettings>
#include <QFileInfo>
#include <QDir>
#include <QFontDatabase>
#include <QMessageBox>
#include <vtkAutoInit.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkObject.h>
#include <vtkOutputWindow.h>
#include <itkOutputWindow.h>
#include "vtkOutputWindow.h"
#include "uStatus.h"
#include "ufunmap.h"
#include "ufunbaseinit.h"
#include "SqlDemo.h"
#include "mainframe.h"
#include "runtimeutil.h"
#include <QMessageBox>
// mitk
#include <QmitkRegisterClasses.h>
#include <mitkRenderingManager.h>
#include <mitkStandaloneDataStorage.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)
VTK_MODULE_INIT(vtkRenderingFreeType)

//查找进程
bool FindProcess(const char * name)
{
    //  Mutex Object
    HANDLE hMutex  =  NULL;
    LPSTR lpszName  = LPSTR(name);

    hMutex  =  ::CreateMutex(NULL, FALSE, lpszName);
    DWORD dwRet  =  ::GetLastError();
    if  (ERROR_ALREADY_EXISTS  ==  dwRet)
    {
        CloseHandle(hMutex);
        return  true;
    }
    return false;
}

int main(int argc, char *argv[])
{
    vtkOutputWindow::SetGlobalWarningDisplay(0);
    QApplication a(argc, argv);
    a.setApplicationName("MitkDemo");
    a.setOrganizationName("None");

    if(FindProcess("MitkDemo.exe"))
    {
        QMessageBox::warning(nullptr,QObject::tr("警告"),QObject::tr("请不要重复运行程序！"));
        return 0;
    }
    if (RuntimeUtil::GetRunning())
    {
        QMessageBox::warning(nullptr, QObject::tr("警告"), QObject::tr("上一次程序运行出现异常结束，请联系供应商！"));
        RuntimeUtil::SetRunning(false);
    }
    QmitkRegisterClasses();
    //
    uFunMap::mFunBaseMap = new QMap<QString, QPointer<uFunBase> >();
    uStatus::mDataStorage = mitk::StandaloneDataStorage::New();

    //
    uStatus::mMainDataDir = QCoreApplication::applicationDirPath() + "/data";
    uStatus::mMainPath = QCoreApplication::applicationDirPath();

    //
    uStatus::mSql=new SqlDemo();

    if(!uStatus::mSql->Initialize())
    {
        return 0;
    }
    uFunBaseInit::getInStance()->f_InitBaseWidgetName();
    MainFrame m;
    m.showMaximized();
    return a.exec();
}
