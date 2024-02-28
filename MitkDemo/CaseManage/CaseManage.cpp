#include "CaseManage.h"
#include "ui_CaseManage.h"
// Qt
#include "QMessageBox"
#include <QTextCodec>
#include <QFileDialog>
#include <QDebug>
#include <QFileDialog>
#include <QPainter>
#include <QDate>
#include <QSize>
#include <QUuid>


// VTK
#include "vtkSmartPointer.h"
#include "vtkWindowToImageFilter.h"
#include "vtkBMPWriter.h"

//

#include "global.h"
#include "systemutil.h"
#include "ufunction.h"
#include "umainfunbase.h"
#include "StdMultiWidget.h"
#include "ustatus.h"
#include "CaseManageControl.h"
#include "SubToolBar.h"
#include "globalsignal.h"
#include <QtConcurrent/QtConcurrent>
#include "progressutil.h"
#include <QJsonDocument>
#include <QJsonObject>
#include <QFuture>
#include <QFutureWatcher>
// mitk
#include <mitkImage.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkProperties.h>
#include "mitkPointSet.h"
#include "QmitkIOUtil.h"

int CaseManage::typeId = qRegisterMetaType<CaseManage*>();

casetable_tuple CaseManage::currentCase;
bool CaseManage::modifyCase = false;
bool CaseManage::isOpen = false;

CaseManage* mPatient;

void f_RefreshPatientFunction()
{
    mPatient->f_refreshTableWidget();
}

void CaseManage::setCheckBox(int pRowIndex, bool isChecked) {
    if (pRowIndex < 0) {
        return;
    }
    if (pRowIndex >= this->ui->tableWidget->rowCount()) {
        return;
    }
    QTableWidgetItem* ckBoxItem =
        static_cast<QTableWidgetItem *>(this->ui->tableWidget->item(pRowIndex, 0));
    if (ckBoxItem == nullptr) {
        ckBoxItem = new QTableWidgetItem();
        this->ui->tableWidget->setItem(pRowIndex, 0, ckBoxItem);
    }
    ckBoxItem->setCheckState(isChecked ? Qt::Checked : Qt::Unchecked);
}

CaseManage::CaseManage(QWidget *parent) :
    uFunBase(parent),
    ui(new Ui::CaseManage)
{
    ui->setupUi(this);
    if(!mSql)
    {
        mSql=new SqlDemo;
    }
    uStatus::pFun_RefreshPatient = &(f_RefreshPatientFunction);

    this->f_refreshTableWidget();

    uStatus::m_uDicomSeries = new THADicomSeries();
    uStatus::mDicomTag = new DicomTag();

    mPatient = this;
    //双击打开病案
    connect(ui->tableWidget,&QTableWidget::itemDoubleClicked,this,&CaseManage::f_openCase);
    connect(ui->tableWidget, &QTableWidget::currentItemChanged, &GlobalSignal::instance(), &GlobalSignal::slot_updateDicom, static_cast<Qt::ConnectionType>(Qt::UniqueConnection | Qt::DirectConnection));
    connect(this,&CaseManage::signalFinishedOpen,this,[&]
    {
        uMainFunBase* mMainFunBase = uFunction::getInStance()->f_GetMain();

        mMainFunBase->f_OpenControl_toobar("SubToolBar");
        dynamic_cast<SubToolBar *>(uFunction::getInStance()->f_GetObjectInstance("SubToolBar*","").data())->f_Refresh();
        mMainFunBase->f_Open_Center("StdMultiWidget");
        QPointer<uFunBase> mStd = uFunction::getInStance()->f_GetObjectInstance("StdMultiWidget*","");
        if (mStd != nullptr)
        {
            dynamic_cast<StdMultiWidget *>(mStd.data())->f_Reset();
        }

    });
}

CaseManage::~CaseManage()
{
    delete ui;
    if(mSql)
    {
        delete mSql;
        mSql=nullptr;
    }
}

QString CaseManage::getName()
{
    return "";
}

QString CaseManage::getID()
{
    return "";
}

bool CaseManage::getClose()
{
    return false;
}

void CaseManage::setCaseOpened(bool b)
{

}

void CaseManage::on_toolButton_save_clicked()
{

}

void CaseManage::f_refreshTableWidget()
{
    QList<casetable_tuple> list_tuple;
    if(gCurrentUser.id!=LOGIN_ID_ADMIN)
    {
        list_tuple=uStatus::mSql->GetCaseTupleByCaseOperator(gCurrentUser.id);
    }
    else
    {
        list_tuple=uStatus::mSql->GetAllCaseTuple();
    }
    f_refreshTableWidget(list_tuple.toVector().toStdVector());
}

void CaseManage::f_refreshTableWidget(const std::vector<casetable_tuple> &pTuples)
{
    ui->tableWidget->horizontalHeader()->setDefaultAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    ui->tableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);
    ui->tableWidget->horizontalHeader()->setVisible(true);
    ui->tableWidget->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidget->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->tableWidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidget->horizontalHeader()->setHighlightSections(false);
    ui->tableWidget->horizontalHeader()->setFixedHeight(40);
    ui->tableWidget->verticalHeader()->setDefaultSectionSize(40);

    ui->tableWidget->clear();
    ui->tableWidget->setRowCount(0);
    int colCount = 8;
    ui->tableWidget->setColumnCount(colCount);
    QTextCodec *codec = QTextCodec::codecForName("utf-8");

    ui->tableWidget->setHorizontalHeaderLabels(QStringList()
        << codec->toUnicode("状态 ")
        << codec->toUnicode("序号")
        << codec->toUnicode("患者姓名 ")
        << codec->toUnicode("患者ID")
        << codec->toUnicode("性别")
        << codec->toUnicode("年龄")
        << codec->toUnicode("医生账户")
        << codec->toUnicode("时间")
        << codec->toUnicode("类型")
       );
    for (size_t i = 0; i < pTuples.size(); i++)
    {
        ui->tableWidget->insertRow(0);
        ui->tableWidget->setItem(0, 0, new QTableWidgetItem(""));
        ui->tableWidget->setItem(0, 1, new QTableWidgetItem(QString::number(pTuples.at(i).case_id)));
        ui->tableWidget->setItem(0, 2, new QTableWidgetItem((pTuples.at(i).patient_name)));
        ui->tableWidget->setItem(0, 3, new QTableWidgetItem((pTuples.at(i).patient_id)));
        ui->tableWidget->setItem(0, 4, new QTableWidgetItem((pTuples.at(i).patient_gender)));
        ui->tableWidget->setItem(0, 5, new QTableWidgetItem(QString::number(pTuples.at(i).patient_age)));
        ui->tableWidget->setItem(0, 6, new QTableWidgetItem((pTuples.at(i).surgery_operator)));
        ui->tableWidget->setItem(0, 7, new QTableWidgetItem((pTuples.at(i).case_create_date.toString("yyyy-MM-dd hh:mm:ss"))));
        setCheckBox(0, 0);
    }
    for (size_t i = 0; i < pTuples.size(); i++)
    {
        for (int j = 0; j < colCount; j++)
        {
            ui->tableWidget->item(i, j)->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
        }
    }
}

void CaseManage::f_appendItems(QStringList pItems)
{
    int columnCount=ui->tableWidget->columnCount();
    if(pItems.size()==0 || pItems.size() % (columnCount) !=0)
        return;
    int index=0;
    int curRow=ui->tableWidget->rowCount()-1;
    for(int i=0;i<pItems.size() ;i++)
    {
        if(i%columnCount==0)
            ui->tableWidget->insertRow(++curRow);
        QTableWidgetItem *item=new QTableWidgetItem(pItems[index]);
        item->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
        ui->tableWidget->setItem(curRow,index%columnCount,item);
        index++;
    }
    ui->tableWidget->selectRow(curRow);
}

void CaseManage::f_setCurrentRow(int pCurRow)
{
    if(pCurRow<0 || pCurRow>=ui->tableWidget->rowCount())
        return;
//    for (int i = 0; i <= ui->tableWidget->rowCount(); i++)
//    {
//        ui->tableWidget->setItem(i,0,new QTableWidgetItem(""));
//    }
    ui->tableWidget->item(pCurRow,0)->setTextAlignment(Qt::AlignCenter);
//    ui->tableWidget->item(pCurRow,0)->setText("*");
    setCheckBox(pCurRow, 1);
    ui->tableWidget->selectRow(pCurRow);
}

void CaseManage::f_setNullSelectRow()
{
    for (int i = 0; i <= ui->tableWidget->rowCount(); i++)
    {
        ui->tableWidget->setItem(i, 0, new QTableWidgetItem(""));
    }
}


void CaseManage::openCase(QString info)
{

}

void CaseManage::closeCase()
{

}


bool CaseManage::f_addCase()
{
    caseinfodialog caseInfo(this);
    // 新增模式
    caseInfo.mIsNew = 1;
    caseInfo.mUUID = "";

    if (QDialog::Rejected == caseInfo.exec())
    {
        return false;
    }
    //刷新列表
    this->f_refreshTableWidget();
    f_setCurrentRow(ui->tableWidget->rowCount());
    mCurrentRowNumbe = ui->tableWidget->rowCount();

    std::vector< std::string >  mPaths;
    QStringList mDcmFiles;

    for (int i = 0; i < uStatus::mDicomSeriesFiles->GetNumberOfTuples(); i++)
    {
        qDebug() << QString::fromStdString( uStatus::mDicomSeriesFiles->GetValue(i) );
        mPaths.push_back(uStatus::mDicomSeriesFiles->GetValue(i));
        mDcmFiles.append(QString(uStatus::mDicomSeriesFiles->GetValue(i)));
    }

    // 加载
//    f_openCaseByID(caseInfo.mCurrentPatientID);
    mCurrentRowNumbe = ui->tableWidget->rowCount() - 1;
    f_setCurrentRow(mCurrentRowNumbe);
    return true;
}


void CaseManage::f_modiCase()
{
    QList<QTableWidgetItem*> list_item = ui->tableWidget->selectedItems();
    if (list_item.size() == 0)
    {
        //System::ER("病例选择页面打开之前未选择病人信息！");
        QMessageBox::information(this, tr("提示"), tr("请先选择要编辑的病案信息！"));
        return;
    }

    int case_id = ui->tableWidget->item(list_item.at(0)->row(), 1)->text().toInt();

    caseinfodialog caseInfo(this);
    // 修改模式
    caseInfo.mIsNew = 0;
    caseInfo.mCurrentPatientID = case_id;
    caseInfo.f_dispCase(case_id);
    caseInfo.SetTitle("修改病例");
    if (QDialog::Rejected == caseInfo.exec())
    {
        return;
    }
    // 处理值
    f_refreshTableWidget();
    //f_setNullSelectRow();
    f_setCurrentRow(mCurrentRowNumbe);
}


void CaseManage::f_openCase()
{
    qDebug()<<"Patient OpenCase start";
    // 判断是否有打开的病例
    if (!uStatus::mCurrentPatientID.isEmpty())
    {
        QMessageBox::warning(this,"提示","请关闭当前病案！");
        return;
    }
    QList<QTableWidgetItem *> items=ui->tableWidget->selectedItems();
    if(items.size()==0)
    {
        QMessageBox::warning(this,"提示","请先选中病案再打开！");
        return;
    }
    int caseId=items.at(CASE_ID)->text().toInt();
    bool isOK{false};
    casetable_tuple tuple=uStatus::mSql->GetCaseTuple(caseId,&isOK);
    if(!isOK)
        return;
    if(tuple.case_dir.isEmpty())
    {
        QMessageBox::warning(this,"提示","请先导入要打开的病案的图像！");
        return;
    }
    gCurrentCase=tuple;
    uStatus::mCurrentPatientID = QString::number(caseId);

    QPointer<uFunBase> b=uFunction::getInStance()->f_GetObjectInstance("StdMultiWidget*","");
    if(b)
    {
        auto w=dynamic_cast<StdMultiWidget *>(b.data());
        if(w)
            uStatus::mMultiWidget=w->GetMultiWidget();
    }
    qDebug() << "Patient::f_openCase 0";
    QList<QTableWidgetItem*> list_item = ui->tableWidget->selectedItems();
    if (list_item.size() == 0)
    {
        //System::ER("病例选择页面打开之前未选择病人信息！");
        QMessageBox::information(this, "提示", "请先选择要打开的病案信息！");
        return;
    }
    mCurrentRowNumbe = ui->tableWidget->selectedItems().at(0)->row();
    f_setCurrentRow(mCurrentRowNumbe);
    f_openCaseByID(caseId);
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance("CaseManageControl*","");
    CaseManageControl*c=dynamic_cast<CaseManageControl*>(mCur.data());
    if (mCur == nullptr || c == nullptr)
        return;
    QMetaObject::invokeMethod(c,"enableButtons",Q_ARG(bool,false));
    qDebug()<<"Patient OpenCase end";
}

void CaseManage::f_closeCase()
{
    qDebug()<<"Patient f_closeCase start";
    qApp->setOverrideCursor(Qt::WaitCursor);
    if (!uStatus::mCurrentPatientID.isEmpty())
    {
        this->f_closeProject();
    }
    qApp->restoreOverrideCursor();
    qDebug()<<"Patient OpenCase end";
}

void CaseManage::f_openCaseByID(int pID)
{
    f_LoadData(pID);
}

void CaseManage::f_LoadData(int )
{
    if (mPlanImageDataMainSeries != nullptr)
    {
        QMessageBox::warning(this,"提示","上次打开的病案未关闭！");
        return;
    }
    // 主序列图像
    QPointer<ProgressUtil> progressUtil = new ProgressUtil;
    progressUtil->Start(tr("正在打开数据"));
    QString dirPath = uStatus::mMainDataDir + "/" + uStatus::mCurrentPatientID ;
    QDir dir(dirPath);
    int step=20;
    
    progressUtil->Step(5);

    if(!QFile(gCurrentCase.case_dir).exists())
    {
        QMessageBox::warning(this,"提示","病案图像丢失，请重新导入！");
        progressUtil->Step(100);
        progressUtil->Finish();
        return;
    }
    mPlanImageDataMainSeries = vtkSmartPointer<Actor>(Actor::New());
    mPlanImageDataMainSeries->f_SetDataStorage(uStatus::mDataStorage);
    mPlanImageDataMainSeries->f_LoadActor();
    progressUtil->Step(step);


    progressUtil->Step(60);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(uStatus::mDataStorage);
    progressUtil->Step(95);

    emit signalFinishedOpen();
    progressUtil->Step(100);
    progressUtil->Finish();
}

void CaseManage::f_UnloadData()
{
    // 主序列图像
    if (mPlanImageDataMainSeries != nullptr)
    {
        mPlanImageDataMainSeries->Delete();
        mPlanImageDataMainSeries=nullptr;
    }

    auto all=uStatus::mDataStorage->GetAll();
    for(auto iter=all->begin();iter!=all->end();iter++)
    {
        uStatus::mDataStorage->Remove(iter->GetPointer());
    }
    uStatus::mCurrentPatientID="";
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void CaseManage::f_LoadProject(int pPatientID)
{

}

void CaseManage::f_closeProject()
{
    f_UnloadData();
    uStatus::mCurrentPatientID="";
}

void CaseManage::f_DelAll()
{
    // 主序列图像
    if (mPlanImageDataMainSeries != nullptr)
    {
        mPlanImageDataMainSeries->Delete();
        mPlanImageDataMainSeries = nullptr;
    }
}

int CaseManage::createDirectory(const std::string &directoryPath)
{
    int dirPathLen = directoryPath.length();
    if (dirPathLen > MAX_PATH_LEN)
    {
        return -1;
    }
    char tmpDirPath[MAX_PATH_LEN] = { 0 };
    for (uint32_t i = 0; i < dirPathLen; ++i)
    {
        tmpDirPath[i] = directoryPath[i];
        if (tmpDirPath[i] == '\\' || tmpDirPath[i] == '/')
        {
            if (ACCESS(tmpDirPath, 0) != 0)
            {
                int32_t ret = MKDIR(tmpDirPath);
                if (ret != 0)
                {
                    return ret;
                }
            }
        }
    }
    return 0;
}

bool CaseManage::f_RemoveCaseByCaseId(int pCaseId)
{
    bool ret=mSql->RemoveCaseByCaseId(QString::number(pCaseId));
    if(!ret)
    {
        return false;
    }
    QList<QTableWidgetItem *> items=ui->tableWidget->selectedItems();
    QString msg="%0删除如下病案：病案序号：%1，患者姓名：%2，患者ID：%3，性别：%4，年龄：%5，医生账户：%6，时间：%7";
    msg=msg.arg(gCurrentUser.id);
    for(int i=1;i<items.size();i++)
    {
        msg=msg.arg(items[i]->text());
    }
    ui->tableWidget->removeRow(ui->tableWidget->currentRow());
    return true;
}

int CaseManage::f_GetSelectedCaseId()
{
    QList<QTableWidgetItem *> items=ui->tableWidget->selectedItems();
    if(items.size()==0)
        return -1;
    return items.at(CASE_ID)->text().toInt();
}

void CaseManage::on_tableWidget_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn)
{
    qDebug() << "CurrentCellChanged:" << currentRow << currentColumn << previousRow << previousColumn;
    if(previousRow<0)
        return;
    setCheckBox(previousRow, 0);
    setCheckBox(currentRow, 1);
    int caseId = ui->tableWidget->item(currentRow,CASE_ID)->text().toInt();
    bool isOK{ false };
    casetable_tuple tuple = uStatus::mSql->GetCaseTuple(caseId, &isOK);
    if (!isOK)
        return;
    QString patientName = tuple.patient_name;
    QString patientId = tuple.patient_id;
    QString text = QString("患者ID:%1  患者姓名：%2").arg(patientId).arg(patientName);
    uStatus::mMain->setPatientInfo(text);
}
