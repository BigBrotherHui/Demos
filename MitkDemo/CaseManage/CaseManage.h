
#ifndef Patient_H
#define Patient_H
#pragma execution_character_set("utf-8")
#include <memory>
#include <QPointer>
#include <QCoreApplication>
#include <QFile>
#include <QDir>
#include <QTimer>
#include <QWidget>


//
#include "SqlDemo.h"
#include "thadicomseries.h"
#include "DicomTag.h"

#define MAX_PATH_LEN 256

#ifdef WIN32
#include <io.h>
#include <direct.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#endif

#ifdef WIN32
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif

#define MAX_PATH_LEN 256

#ifdef WIN32
#include <io.h>
#include <direct.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#endif

#ifdef WIN32
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif



//
#include "ufunbase.h"
//
#include "caseinfodialog.h"
#include <QHeaderView>
#include "Actor.h"
class QTableWidgetItem;
class CaseInfoDialog;
class SqlDemo;
namespace Ui {
class CaseManage;
}

class CaseManage : public uFunBase
{
    Q_OBJECT

public:
    enum CASE_
    {
        CASE_STATE=0,
        CASE_ID,
        CASE_PATIENTNAME,
        CASE_PATIENTID,
        CASE_PATIENTGENDER,
        CASE_PATIENTAGE,
        CASE_DOCTORNAME,
    };
    static int typeId;
    Q_INVOKABLE explicit CaseManage(QWidget *parent = nullptr);
    ~CaseManage();

    //
    static const std::string PROJECT_PATH;
    static casetable_tuple currentCase;
    static QString currentMITKFileName;
    static bool modifyCase;

    void setCheckBox(int pRowIndex, bool isChecked);
    QString getName();
    QString getID();
    bool getClose();
    static bool isOpen;
    void setCaseOpened(bool b);
    //
    void f_refreshTableWidget();
    void f_refreshTableWidget(const std::vector<casetable_tuple> &pTuples);
    void f_appendItems(QStringList pItems);
    void f_setCurrentRow(int pCurRow);
    void f_setNullSelectRow();

    bool f_addCase();
    void f_modiCase();
    void f_openCase();
    void f_closeCase();
    bool f_RemoveCaseByCaseId(int pCaseId);
    int  f_GetSelectedCaseId();
  Q_SIGNALS:
      void signalFinishedOpen();

private:

    QMetaObject::Connection eidtConnect;
    void openCase(QString info);
    void closeCase();

    itk::WeakPointer<mitk::DataNode> outputNode;

    QAction* closeProjectAction;

    // 建目录
    int createDirectory(const std::string &directoryPath);

    int mCurrentRowNumbe = 0;

    vtkSmartPointer<Actor> mPlanImageDataMainSeries{nullptr};

private Q_SLOTS:
    void on_toolButton_save_clicked();
    void on_tableWidget_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn);

protected:

    void f_openCaseByID(int pID);
    void f_LoadData(int pID);
    void f_UnloadData();
    void f_closeProject();
    void f_DelAll();

    void f_LoadProject(int pPatientID);

private:
    Ui::CaseManage *ui;
    SqlDemo *mSql{nullptr};
    int mLastSelectedRow{-1};
};

#endif
