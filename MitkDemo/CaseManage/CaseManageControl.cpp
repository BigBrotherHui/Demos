#include <QDebug>
#include <QVBoxLayout>
#include "QMessageBox"
#include <QtConcurrent/QtConcurrent>
#include <QFileDialog>
#include <QJsonDocument>
#include <QJsonObject>
#include "ustatus.h"
#include "umainfunbase.h"
#include "CaseManageControl.h"
#include "ui_CaseManageControl.h"
#include "CaseManage.h"
#include "StdMultiWidget.h"
#include "vtkInteractorStyleImage.h"
#include "systemutil.h"
#include <QTableWidgetItem>
#include "global.h"
#include "caseimagedialog.h"
#include "QMessageBox.h"
#include <vtkRenderer.h>
#include "globalsignal.h"

vtkSmartPointer<vtkDICOMImageReader>CaseManageControl::RenderInfo::reader=vtkSmartPointer<vtkDICOMImageReader>::New();
int CaseManageControl::typeId = qRegisterMetaType<CaseManageControl*>();
CaseManageControl::CaseManageControl(QWidget *parent) :
    uFunBase(parent),
    ui(new Ui::CaseManageControl)
{
    ui->setupUi(this);
    mRenderTimer=new QTimer(this);
    mRenderTimer->setInterval(50);
    if(!mImageViewer)
    {
        mImageViewer=vtkSmartPointer<vtkImageViewer2>::New();
        mImageViewer->GetRenderer()->GetActiveCamera()->SetParallelProjection(1);
        ui->widget_dicomWidget->interactor()->SetInteractorStyle(vtkInteractorStyleImage::New());
        ui->widget_dicomWidget->renderWindow()->AddRenderer(mImageViewer->GetRenderer());
        mImageViewer->GetRenderWindow()->SetOffScreenRendering(1);
    }
    connect(mRenderTimer,&QTimer::timeout,this,[&]
    {
        if(mRenderInfoVec.size()==0)
        {
            mRenderTimer->stop();
            return;
        }
        qApp->setOverrideCursor(Qt::WaitCursor);
        RenderInfo info=mRenderInfoVec.top();
        while(mRenderInfoVec.size()>0)
            mRenderInfoVec.pop();

        if(info.imageDir.isEmpty())
        {
            mImageViewer->SetInputData(vtkImageData::New());
            ui->verticalSlider_slice->blockSignals(1);
            ui->verticalSlider_slice->setValue(info.slice);
            ui->verticalSlider_slice->blockSignals(0);
            ui->widget_dicomWidget->renderWindow()->Render();
            qApp->restoreOverrideCursor();
            return;
        }            

        info.reader->SetDirectoryName(info.imageDir.toStdString().c_str());
        info.reader->Update();
        mImageViewer->SetInputData(info.reader->GetOutput());
        ui->verticalSlider_slice->setRange(mImageViewer->GetSliceMin(),mImageViewer->GetSliceMax());
        mImageViewer->GetRenderer()->ResetCamera();
        ui->verticalSlider_slice->blockSignals(1);
        ui->verticalSlider_slice->setValue(info.slice);
        ui->verticalSlider_slice->blockSignals(0);
        mImageViewer->SetSlice(info.slice);
        double scale=mImageViewer->GetRenderer()->GetActiveCamera()->GetParallelScale();
        mImageViewer->GetRenderer()->ResetCamera();
        mImageViewer->GetRenderer()->GetActiveCamera()->SetParallelScale(scale);
        ui->widget_dicomWidget->renderWindow()->Render();

        qApp->restoreOverrideCursor();
    });
    connect(&GlobalSignal::instance(), &GlobalSignal::signal_updateDicom, this, &CaseManageControl::slot_updateDicom, static_cast<Qt::ConnectionType>(Qt::DirectConnection | Qt::UniqueConnection));
    mEnableButtons.append(ui->cb_del);
    mEnableButtons.append(ui->pushButton_image);
    mEnableButtons.append(ui->cb_new);
    mEnableButtons.append(ui->cb_modi);
    mEnableButtons.append(ui->cb_open);
    mEnableButtons.append(ui->cb_import);
    mEnableButtons.append(ui->cb_export);
    ui->cb_close->setEnabled(0);
    ui->checkBox_previewDicom->SetCheckedBkColor("#1167CF");
    ui->checkBox_previewDicom->SetCheckedHandleColor(Qt::white);
    ui->checkBox_previewDicom->SetUncheckedBkColor("#131720");
    ui->checkBox_previewDicom->SetUncheckedHandleColor(Qt::white);
}

CaseManageControl::~CaseManageControl()
{
    delete ui;
}

void CaseManageControl::ShowDicom()
{
    RenderInfo renderInfo;
    renderInfo.imageDir=mCurrentImageDir;
    renderInfo.slice=0;
    mRenderInfoVec.push(renderInfo);
    if(!mRenderTimer->isActive())
        mRenderTimer->start();
}

void CaseManageControl::slot_updateDicom(QTableWidget* w, QTableWidgetItem* cur, QTableWidgetItem* prev)
{
    if (!cur || !w || (w->selectedItems().size() == 0 && !prev))
        return;
    int curRow = cur->row();
    if (curRow == mLastSelectedRow)
        return;
    mLastSelectedRow = curRow;
    int caseId = w->item(curRow, CaseManage::CASE_ID)->text().toInt();
    bool isOK{ false };
    casetable_tuple ca = uStatus::mSql->GetCaseTuple(caseId, &isOK);
    if (!isOK)
        return;
    QString imageDir = ca.case_dir;
    if (!QDir(imageDir).exists())
    {
        qDebug() << QString("dicom路径：%1不存在！").arg(imageDir);
        return;
    }
    mCurrentImageDir = imageDir;
    if (!ui->checkBox_previewDicom->isChecked())
        return;
    ShowDicom();
}

void CaseManageControl::enableButtons(bool pIsEnable)
{
    for(int i=0;i<mEnableButtons.size();i++)
    {
        mEnableButtons.at(i)->setEnabled(pIsEnable);
    }
    if(pIsEnable)
        ui->cb_close->setEnabled(false);
    else
        ui->cb_close->setEnabled(1);
}

void CaseManageControl::on_cb_new_clicked()
{
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance("CaseManage*","");
    if (mCur != nullptr)
    {
        bool ret=dynamic_cast<CaseManage*>(mCur.data())->f_addCase();
        if(ret && ui->checkBox_previewDicom->isChecked())
            on_checkBox_previewDicom_clicked();
        if(ret)
        {
            if(QMessageBox::Yes == QMessageBox::question(nullptr,"提示","病案创建完成，是否为病案添加图像？", QMessageBox::Yes | QMessageBox::No))
                on_pushButton_image_clicked();
        }
    }
}

void CaseManageControl::on_cb_modi_clicked()
{
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance("CaseManage*","");
     if (mCur != nullptr)
     {
         dynamic_cast<CaseManage*>(mCur.data())->f_modiCase();
     }
}

void CaseManageControl::on_cb_open_clicked()
{
    if(!uStatus::mCurrentPatientID.isEmpty())
    {
        QMessageBox::warning(this,"提示","请先选中病案再打开！");
        return;
    }
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance("CaseManage*","");
    CaseManage*patient=dynamic_cast<CaseManage*>(mCur.data());
    if (mCur == nullptr || patient == nullptr)
        return;
    patient->f_openCase();
}

void CaseManageControl::on_pushButton_search_clicked()
{
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance("CaseManage*","");
    CaseManage*patient=dynamic_cast<CaseManage*>(mCur.data());
    if (mCur == nullptr || patient == nullptr)
        return;
    QString keyword=ui->lineEdit_search->text();
    if(keyword.isEmpty())
    {
        patient->f_refreshTableWidget();
        return;
    }
    std::vector<casetable_tuple> ret;
    if(gCurrentUser.id==LOGIN_ID_ADMIN)
    {
        ret=uStatus::mSql->SearchCaseByKeyword(keyword);
    }
    else
    {
        ret=uStatus::mSql->SearchCaseByKeywordAndCaseOperator(keyword,gCurrentUser.id);
    }
    patient->f_refreshTableWidget(ret);
}

void CaseManageControl::on_lineEdit_search_returnPressed()
{
    on_pushButton_search_clicked();
}

void CaseManageControl::on_cb_del_clicked()
{
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance("CaseManage*","");
    CaseManage*patient=dynamic_cast<CaseManage*>(mCur.data());
    if (mCur == nullptr || patient == nullptr)
        return;
    int caseId=patient->f_GetSelectedCaseId();
    if(caseId==-1)
    {
        QMessageBox::warning(this,"提示","请先选中待删除病案！");
        return;
    }
    bool isOk{false};
    casetable_tuple tuple=uStatus::mSql->GetCaseTuple(caseId,&isOk);
    if (QMessageBox::No == QMessageBox::question(nullptr, "提示","确认要删除选中的病案？",QMessageBox::Yes | QMessageBox::No))
    {
        return;
    }
    bool ret=patient->f_RemoveCaseByCaseId(caseId);
    if(!ret)
    {
        QMessageBox::warning(this,"提示","删除失败！");
        return;
    }
    QString fullPath=uStatus::mMainDataDir+"/"+QString::number(caseId);
    ret=SystemUtil::DeleteDir(fullPath);
    qDebug()<<fullPath;
    if(ret==false)
    {
        QMessageBox::warning(this,"提示","删除失败！");
    }
    else
    {
        QMessageBox::warning(this,"提示","删除成功！");
    }
}

void CaseManageControl::on_verticalSlider_slice_valueChanged(int value)
{
    if(!mImageViewer->GetInput())
        return;
    mImageViewer->SetSlice(value);
    double scale=mImageViewer->GetRenderer()->GetActiveCamera()->GetParallelScale();
    mImageViewer->GetRenderer()->ResetCamera();
    mImageViewer->GetRenderer()->GetActiveCamera()->SetParallelScale(scale);
    ui->widget_dicomWidget->renderWindow()->Render();
}

void CaseManageControl::on_checkBox_previewDicom_clicked()
{
    bool isChecked=ui->checkBox_previewDicom->isChecked();
    RenderInfo info;
    if(isChecked)
    {
        info.imageDir=mCurrentImageDir;
        info.scale=1.0;
    }
    info.slice=ui->verticalSlider_slice->value();
    mRenderInfoVec.push(info);
    if(!mRenderTimer->isActive())
    {
        mRenderTimer->start();
    }
}

void CaseManageControl::on_cb_close_clicked()
{
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance("CaseManage*","");
    CaseManage*patient=dynamic_cast<CaseManage*>(mCur.data());
    if(mCur==nullptr || patient==nullptr)
    {
        return;
    }
    uMainFunBase* mMainFunBase = uFunction::getInStance()->f_GetMain();
    mMainFunBase->f_CloseControl_right("uBlank");

    mMainFunBase->f_CloseControl_Center("StdMultiWidget");

    patient->f_closeCase();

    mMainFunBase->f_OpenControl_toobar("MainToolBar");
    enableButtons(1);
}

void CaseManageControl::on_pushButton_image_clicked()
{
    QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance("CaseManage*","");
    CaseManage*patient=dynamic_cast<CaseManage*>(mCur.data());
    if(mCur==nullptr || patient==nullptr)
    {
        return;
    }
    int id=patient->f_GetSelectedCaseId();
    if(id==-1)
    {
        QMessageBox::warning(this,"提示","请先选中病案！");
        return;
    }
    CaseImageDialog d;
    d.SetCurrentCaseId(id);
    if(QDialog::Accepted==d.exec())
    {
        bool isChecked=ui->checkBox_previewDicom->isChecked();
        if(isChecked)
        {
            RenderInfo renderInfo;
            QPointer<uFunBase> mCur = uFunction::getInStance()->f_GetObjectInstance("CaseManage*","");
            CaseManage*patient=dynamic_cast<CaseManage*>(mCur.data());
            if (mCur == nullptr || patient == nullptr)
                return;
            int caseId=patient->f_GetSelectedCaseId();
            if(caseId==-1)
                return;
            bool ok{false};
            casetable_tuple t=uStatus::mSql->GetCaseTuple(caseId,&ok);
            if(!ok)
                return;
            mCurrentImageDir=t.case_dir;
            renderInfo.imageDir=mCurrentImageDir;
            renderInfo.slice=0;
            mRenderInfoVec.push(renderInfo);
            if(!mRenderTimer->isActive())
                mRenderTimer->start();
        }
    }
}
