#pragma execution_character_set("utf-8")
#include "ustatus.h"
#include <QSettings>
#include <vtkDICOMReader.h>

#include "QmitkStdMultiWidget.h"
#define tr QObject::tr
MainWindow *uStatus::mMain=nullptr;
QString uStatus::mVer = "V1";
QString uStatus::mKernelVer = "V1.0.0.0.88";
QString uStatus::mImageObjectPathName;
QMessageLogContext uStatus::context;
//Main
QString uStatus::mMainDataDir="";
QString uStatus::mLogDir="log";
QString uStatus::mBackupDir="backup";
QString uStatus::mDicomTmpDir="dicomTmpDir";
QString uStatus::mMainConfigureFileName = "syscfg.ini";
QString uStatus::mMainWindow="MainWindow";
QString uStatus::mMainToolbarWidget="uBlank";
QString uStatus::mMainCenterWidget="uBlank";
QString uStatus::mMainControlWidgetRight="uBlank";
QString uStatus::mMainControlWidgetLeft="uBlank";
QString uStatus::mMainFootbarWidget="NULL";
//
QString uStatus::mMainPath = "";
//QString uStatus::mCurClassName = "";
//QString uStatus::mCurClassID = "";
QString uStatus::mCurrentPatientID;
QString uStatus::mDicomSubDir = "/dcm";
QString uStatus::mDicomOrgDir = "/orgdcm";
QString uStatus::mBMPSubDir = "/bmp";
QString uStatus::mParameterSubDir = "/parameter";
QString uStatus::mExternalGenDir = "";


mitk::DataStorage::Pointer uStatus::mDataStorage = nullptr;

QmitkStdMultiWidget *uStatus::mMultiWidget = nullptr;

SqlDemo *uStatus::mSql= nullptr;
//
int uStatus::mUserLevel = 0;
QString uStatus::mUserCode = "";
THADicomSeries* uStatus::m_uDicomSeries = nullptr;
vtkStringArray* uStatus::mDicomSeriesFiles = nullptr;
vtkDICOMReader* uStatus::mDICOMReader = nullptr;
DicomTag* uStatus::mDicomTag = nullptr; 
int uStatus::mExtent[6] = {0,0,0,0,0,0};
double uStatus::mOrigin[3] = {0,0,0};
double uStatus::mSpacing[3] = {0,0,0};
double uStatus::mCenter[3] = {0,0,0};
int uStatus::Dims[3] = {0,0,0};
int uStatus::mWW=100;
int uStatus::mWL=100;
int uStatus::mDefaultWW=100;
int uStatus::mDefaultWL=100;

//MITK
QAction* uStatus::gActionCloseProject = nullptr;
QToolBar* uStatus::gMainToolBar = nullptr;

void(*uStatus::pFun_AIViewSetCurrent)() = nullptr;
//void(*uStatus::pFun_AIProsthesisIVAngel)() = nullptr;
void(*uStatus::pFun_Navigator)(double pX, double pY, double pZ)=nullptr;
void(*uStatus::pFun_SetProsthesisFixStatus)(int pStatus);
void(*uStatus::pFun_RefreshPatient)();
void(*uStatus::pFun_UpdateHash)(int case_ID, QString key, QVariant value);

