
#ifndef caseinfodialog_H
#define caseinfodialog_H
#pragma execution_character_set("utf-8")
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


#include <QDialog>
#include <QPointer>
#include "vtkSmartPointer.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkImageViewer2.h"
#include "vtkCamera.h"
#include "vtkStringArray.h"
#include "QVTKOpenGLNativeWidget.h"

//
#include "SqlDemo.h"

#include "EmptyImageData.h"

//class CaseSelectionPage;
class CaseManage;

namespace Ui {
class caseinfodialog;
}

class caseinfodialog : public QDialog
{
    Q_OBJECT

public:
    explicit caseinfodialog(QWidget *parent = nullptr);
    ~caseinfodialog();
    void SetInfo(casetable_tuple info);
    void SetEditable();
    casetable_tuple GetInfo();
    int f_createNewCase(casetable_tuple info);
    void f_onCaseupData(casetable_tuple modifiedCase);
    void f_updateCase(casetable_tuple info);
    void f_dispCase(int pPatientID);
    void SetTitle(QString strtitle);
    void keyPressEvent(QKeyEvent *event);
    int caseinfodialog::createDirectory(const std::string &directoryPath);
    CaseManage* caseSelectionPage;

    int mIsNew = 0;
    QString mUUID = "";
    int mCurrentPatientID = 0;


private slots:
    void on_pushButton_done_clicked();

    void on_pushButton_exit_clicked();

    void onCaseupData();

    void on_pushButton_close_clicked();

signals:
    void Caseupdata();

    friend class CaseSelectionPage;
private:
    SqlDemo sql;
    Ui::caseinfodialog *ui;
};

#endif // caseinfodialog_H
