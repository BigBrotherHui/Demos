#pragma execution_character_set("utf-8")
#include "caseimagedialog.h"
#include "ui_caseimagedialog.h"

#include <QFileDialog>
#include "QMessageBox.h"
#include <QDebug>
#include "ustatus.h"
#include "systemutil.h"
#include "thadicomseries.h"
#include "DicomTag.h"
#include <vtkCamera.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include "SqlDemo.h"

CaseImageDialog::CaseImageDialog(QWidget* parent):
    QDialog(parent),
    ui(new Ui::CaseImageDialog)
{
    qDebug() << "CaseImageDialog::CaseImageDialog";
    ui->setupUi(this);
    setWindowFlag(Qt::FramelessWindowHint);
    setWindowFlag(Qt::SubWindow);
    mImageView=vtkSmartPointer<vtkImageViewer2>::New();
    mImageView->GetRenderer()->GetActiveCamera()->SetParallelProjection(1);
    ui->widget_3d->interactor()->SetInteractorStyle(vtkInteractorStyleImage::New());
    ui->widget_3d->renderWindow()->AddRenderer(mImageView->GetRenderer());
    mImageView->GetRenderWindow()->SetOffScreenRendering(1);
    dialog = QPointer<SeriesSelectDialog>(new SeriesSelectDialog(this));
    uStatus::mDicomSeriesFiles = vtkStringArray::New();
    ui->verticalSlider_slice->setEnabled(0);
    qDebug() << "CaseImageDialog::CaseImageDialog finished";
}

CaseImageDialog::~CaseImageDialog()
{

    delete ui;
    if(dialog)
    {
        delete dialog;
        dialog=nullptr;
    }
    qDebug() << "CaseImageDialog::~CaseImageDialog finished";
}

void CaseImageDialog::on_pushButton_image_clicked()
{
    QString filePath = QFileDialog::getExistingDirectory(this, tr("请选择图像路径..."), "D:/");
    if (filePath.isEmpty())
    {
      return;
    }
    if(!QDir(filePath).exists())
    {
        QMessageBox::warning(this,"提示","导入路径不存在！");
        return;
    }
    if(!SystemUtil::CheckDir(uStatus::mMainDataDir+"/"+uStatus::mDicomTmpDir,1))
    {
        QMessageBox::warning(this,"提示",QString("导入失败！%1下没有足够的空间！").arg(uStatus::mMainDataDir));
        return;
    }
    SystemUtil::DeleteDir(uStatus::mMainDataDir+"/"+uStatus::mDicomTmpDir);
    SystemUtil::CopyDir(filePath,uStatus::mMainDataDir+"/"+uStatus::mDicomTmpDir);
    dialog->loadSeries(uStatus::mMainDataDir+"/"+uStatus::mDicomTmpDir);
    if (dialog->exec() == QDialog::Accepted)
    {
        uStatus::mDicomSeriesFiles->DeepCopy(dialog->GetFileNamesForSeries());
        if(!uStatus::m_uDicomSeries)
        {
            uStatus::m_uDicomSeries=new THADicomSeries;
        }
        uStatus::m_uDicomSeries->f_GetDicomSeries(uStatus::mDicomSeriesFiles);
        if(!uStatus::mDicomTag)
        {
            uStatus::mDicomTag=new DicomTag;
        }
        uStatus::mDicomTag->f_SetSystemData();
        mImageView->SetInputData(uStatus::m_uDicomSeries->f_GetImageData());
        mImageView->GetRenderer()->GetActiveCamera()->SetViewUp(0, -1, 0);
        mImageView->GetRenderer()->GetActiveCamera()->SetPosition(0, 0, -1);
        mImageView->GetRenderer()->ResetCamera();

        this->ui->verticalSlider_slice->setMaximum(uStatus::Dims[2]);
        this->ui->verticalSlider_slice->setMinimum(1);

        ui->widget_3d->renderWindow()->Render();
        dialog->deleteLater();
        ui->verticalSlider_slice->setEnabled(1);
    }
}

void CaseImageDialog::on_pushButton_close_clicked() { this->reject(); }

void CaseImageDialog::on_verticalSlider_slice_valueChanged(int pValue)
{
    mImageView->SetSlice(pValue);
    double scale=mImageView->GetRenderer()->GetActiveCamera()->GetParallelScale();
    mImageView->GetRenderer()->ResetCamera();
    mImageView->GetRenderer()->GetActiveCamera()->SetParallelScale(scale);
    ui->widget_3d->renderWindow()->Render();
}

void CaseImageDialog::on_pushButton_import_clicked()
{
    qDebug() << "CaseImageDialog : import button down... ";
    if(mCaseId==-1)
    {
        QMessageBox::warning(this,"提示","未设置病案号！");
        return;
    }
    int imageNum=uStatus::mDicomSeriesFiles->GetNumberOfTuples();
    if(imageNum==0)
    {
        QMessageBox::warning(this,"提示","请先选择图像！");
        return;
    }
    if(QMessageBox::No == QMessageBox::question(this,"提示","确定导入该组图像吗？导入后原图像会被清空！", QMessageBox::Yes | QMessageBox::No))
        return;
    QString mTagFilePath = uStatus::mMainDataDir + "/" + QString::number(mCaseId) + uStatus::mDicomOrgDir;
    SystemUtil::DeleteDir(mTagFilePath);
    for (int i = 0; i < imageNum; i++)
    {
        QFileInfo info = QFileInfo(uStatus::mDicomSeriesFiles->GetValue(i).c_str());
        QString mTagFileName = mTagFilePath + "/" + info.fileName();
        if (!QFile::exists(mTagFilePath))
        {
            QDir().mkpath(mTagFilePath);
        }
        QFile::copy(uStatus::mDicomSeriesFiles->GetValue(i).c_str(), mTagFileName);
    }
    uStatus::mSql->UpdateCaseTableTuple(mCaseId,"case_directory",mTagFilePath);
    SystemUtil::DeleteDir(uStatus::mMainDataDir+"/"+uStatus::mDicomTmpDir);
    QMessageBox::warning(this,"提示","导入成功！");
    this->accept();
}
