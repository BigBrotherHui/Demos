#pragma execution_character_set("utf-8")
#include "ProgressWidget.h"
//#include "ui_ProgressWidget.h"
#include "caseinfodialog.h"
#include "ui_caseinfodialog.h"

#include "thadicomseries.h"
#include "DicomTag.h"
#include <QFileDialog>
#include <QFileInfo>
#include "QMessageBox"
#include <QRegExp>
#include <QString>
#include <QRegExpValidator>
#include <QUuid>
#include <QDesktopWidget>
#include <QKeyEvent>
#include <QProgressBar>

#include "vtkRenderer.h"
#include <vtkInteractorStyleImage.h>
#include "CaseManage.h"
#include "systemutil.h"
#include "global.h"
#include <QListView>
QString getChinese(QString str)
{
  //    QString str;
  QString chineseStr;
  int nCount = str.count();
  for (int i = 0; i < nCount; i++)
  {
    QChar cha = str.at(i);
    ushort uni = cha.unicode();
    if (uni >= 0x4E00 && uni <= 0x9FA5)
    {
      //这个字符是中文
      qDebug() << uni;
      chineseStr.append(uni);
      //            qDebug("%s",uni);
    }
  }
  return chineseStr; //最后返回的这个字符串就是中文字符串
}

caseinfodialog::caseinfodialog(QWidget *parent) : QDialog(parent), ui(new Ui::caseinfodialog)
{
  ui->setupUi(this);
  this->setWindowFlag(Qt::FramelessWindowHint);
  ui->comboBox_patientGender->setView(new QListView(ui->comboBox_patientGender));
  ui->comboBox_userid->setView(new QListView(ui->comboBox_userid));
  ui->widget->closeByOthers();
  ui->widget->setTitle("新建病例");

  ui->lineEdit_patientName->setMaxLength(18);
  ui->lineEdit_patientID->setMaxLength(16);
  if(gCurrentUser.id==LOGIN_ID_ADMIN)
  {
      std::vector<user_info> infos;
      bool ret=uStatus::mSql->GetUserByEnable(infos,Enable_Able);
      if(ret)
      {
          for(size_t i=0;i<infos.size();i++)
          {
              ui->comboBox_userid->insertItem(0,infos.at(i).id);
          }
      }
      ui->comboBox_userid->setEnabled(1);
  }
  else
  {
      ui->comboBox_userid->insertItem(0,gCurrentUser.id);
      ui->comboBox_userid->setEnabled(0);
  }

  //
  QRegExp regx2("[a-zA-Z0-9]+$");
  QValidator *validator_patientID = new QRegExpValidator(regx2, ui->lineEdit_patientID);
  ui->lineEdit_patientID->setValidator(validator_patientID);

  connect(ui->widget, &TitleBar::signal_closeWindow, this, [&]() {
    if (mIsNew == 1)
    {
        if (QMessageBox::Yes == QMessageBox::question(this,"提示","确定要取消新增病案吗？", QMessageBox::Yes | QMessageBox::No))
      {
        close();
      }
    }
    else
    {
      if (QMessageBox::Yes == QMessageBox::question(this,"提示","确定要取消修改病案吗？", QMessageBox::Yes | QMessageBox::No))
      {
        close();
      }
      //close();
    }
  });
  //修改数据需要重新获取数据
  if (mIsNew == 0)
  {
    //显示当前选中数据
  }
  else
  {
    //新增数据
  }
}

caseinfodialog::~caseinfodialog()
{
  delete ui;
}

void caseinfodialog::SetInfo(casetable_tuple info)
{
  ui->comboBox_patientGender->setCurrentText(info.patient_gender);
  // ui->lineEdit_caseID->setText(QString::number(info.case_id));
  ui->lineEdit_patientName->setText(info.patient_name);
  ui->lineEdit_patientID->setText(info.patient_id);
  ui->lineEdit_patientAge->setValue(info.patient_age);
  ui->comboBox_userid->setCurrentText(info.surgery_operator);
}

void caseinfodialog::SetTitle(QString strtitle)
{
  ui->widget->setTitle(strtitle);
}

void caseinfodialog::keyPressEvent(QKeyEvent *event)
{
  int nkey = event->key();
  switch (event->key())
  {
    case Qt::Key_Escape: // 按下的为Esc键
      break;             // 不做反应直接退出
    case Qt::Key_Return: // 按下的为Enter键
    {
      //on_pushButton_browse_clicked();
    }
    break; // 不做反应直接退出
    default:
      QDialog::keyPressEvent(event);
  }
}

void caseinfodialog::on_pushButton_exit_clicked()
{
  this->close();
}

void caseinfodialog::on_pushButton_done_clicked()
{
  // qDebug() << "新建或修改按钮被按下...";
  casetable_tuple modifiedCase = CaseManage::currentCase;

  if (ui->lineEdit_patientID->text().isEmpty())
  {
    qDebug() << "病例信息对话框未输入患者ID";
    QMessageBox::information(this, "提示", "请输入患者ID");
    return;
  }

  if (ui->lineEdit_patientName->text().isEmpty())
  {
    qDebug() << QString::fromLocal8Bit("病例信息对话框未输入患者姓名 ");
    QMessageBox::information(this, "提示 ", "请输入患者姓名");
    return;
  }
  modifiedCase.patient_name = ui->lineEdit_patientName->text();
  modifiedCase.patient_name.trimmed();
  modifiedCase.patient_name.replace(" ", "");

  modifiedCase.patient_id = ui->lineEdit_patientID->text();

  if (ui->lineEdit_patientAge->text().isEmpty())
  {
    qDebug() << "病例信息对话框未输入患者年龄 ";
    QMessageBox::information(this, "提示 ", "请输入患者年龄");
    return;
  }
  modifiedCase.patient_age = ui->lineEdit_patientAge->value();
  if (modifiedCase.patient_age > 149)
  {
    QMessageBox::information(this, "提示 ", "输入患者年龄不能大于150");
    return;
  }

  if (ui->comboBox_userid->currentText().isEmpty())
  {
    qDebug() << "病例信息对话框未输入医生姓名 ";
    QMessageBox::information(this, "提示 ", "请输入医生姓名");
    return;
  }
  modifiedCase.surgery_operator = ui->comboBox_userid->currentText();
  modifiedCase.surgery_operator.trimmed();
  modifiedCase.surgery_operator.replace(" ", "");

  if (ui->comboBox_patientGender->currentText().isEmpty())
  {
    qDebug() << "病例信息对话框未输入患者性别 ";
    QMessageBox::information(this, "提示 ", "请输入患者性别");
    return;
  }

  // 判断新增还是修改
  if (mIsNew == 0)
  {
    // 修改 修改部分信息，不处理dicom文件等信息，只是数据库中的信息的更新
    modifiedCase.case_id = mCurrentPatientID;
    modifiedCase.patient_gender = ui->comboBox_patientGender->currentText();
    f_onCaseupData(modifiedCase);
  }
  else
  {
    modifiedCase.patient_gender = ui->comboBox_patientGender->currentText();
    //新增 保持数据以后，复制文件到指定位置
    modifiedCase.case_id = this->sql.GenerateUniqueCaseID();
    modifiedCase.case_create_date = QDateTime::currentDateTime();
    modifiedCase.case_operator=ui->comboBox_userid->currentText();
    QString mTagFilePath = uStatus::mMainDataDir + "/" + QString::number(modifiedCase.case_id) + uStatus::mDicomOrgDir;
    modifiedCase.case_dir.clear();
    int ret = this->f_createNewCase(modifiedCase);
    mCurrentPatientID = modifiedCase.case_id;
    // 复制文件到制定目录
    if(!uStatus::mDicomSeriesFiles)
    {
        uStatus::mDicomSeriesFiles=vtkStringArray::New();
    }
    int imageNum=uStatus::mDicomSeriesFiles->GetNumberOfTuples();
    if(imageNum!=0)
    {
        QDir tagDir(mTagFilePath);
        if(tagDir.entryList(QDir::Files | QDir::NoDotAndDotDot | QDir::Dirs).size()>0)
            SystemUtil::DeleteDir(mTagFilePath);
        for (int i = 0; i < uStatus::mDicomSeriesFiles->GetNumberOfTuples(); i++)
        {
            QFileInfo info = QFileInfo(uStatus::mDicomSeriesFiles->GetValue(i).c_str());
            QString mTagFileName = mTagFilePath + "/" + info.fileName();
            qDebug() << "Copy Files : " << mTagFileName;
            if (!QFile::exists(mTagFilePath))
            {
                QDir().mkpath(mTagFilePath);
            }
            QFile::copy(uStatus::mDicomSeriesFiles->GetValue(i).c_str(), mTagFileName);
        }
        SystemUtil::DeleteDir(uStatus::mMainDataDir+"/"+uStatus::mDicomTmpDir);
    }
  }
  this->accept();
}

void caseinfodialog::onCaseupData() {}

casetable_tuple caseinfodialog::GetInfo()
{
  casetable_tuple tp;
  tp.patient_name = ui->lineEdit_patientName->text();
  tp.patient_id = ui->lineEdit_patientID->text();
  tp.patient_gender = ui->comboBox_patientGender->currentText();
  tp.patient_age = ui->lineEdit_patientAge->value();
  tp.surgery_operator = ui->comboBox_userid->currentText();
  tp.case_operator=gCurrentUser.id;
  return tp;
}

void caseinfodialog::SetEditable()
{
  ui->comboBox_patientGender->setCurrentText("男 ");
  ui->lineEdit_patientName->setText(" ");
  ui->lineEdit_patientID->setText(" ");
  ui->lineEdit_patientAge->setValue(0);
  ui->comboBox_userid->currentText();
}

int caseinfodialog::createDirectory(const std::string &directoryPath)
{
  int dirPathLen = directoryPath.length();
  if (dirPathLen > MAX_PATH_LEN)
  {
    return -1;
  }
  char tmpDirPath[MAX_PATH_LEN] = {0};
  for (int i = 0; i < dirPathLen; ++i)
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

// 保存
int caseinfodialog::f_createNewCase(casetable_tuple info)
{
  if (this->sql.InsertCaseTableTuple(info))
  {
    return -1;
  }

  return 0;
}

void caseinfodialog::f_onCaseupData(casetable_tuple modifiedCase)
{
  this->f_updateCase(modifiedCase);
}

void caseinfodialog::f_updateCase(casetable_tuple info)
{
  if (this->sql.UpdateCaseTableTuple(info))
  {
    qDebug() << "Fail to update case.";
  }
}

void caseinfodialog::f_dispCase(int pPatientID)
{
  casetable_tuple currentCase = CaseManage::currentCase;
  currentCase = this->sql.GetCaseTuple(pPatientID);

  ui->lineEdit_patientName->setText(currentCase.patient_name);

  ui->lineEdit_patientID->setText(currentCase.patient_id);

  if(currentCase.patient_age<0)
      currentCase.patient_age=0;
  ui->lineEdit_patientAge->setValue(currentCase.patient_age);

  ui->comboBox_userid->setCurrentText(currentCase.surgery_operator);

  ui->comboBox_patientGender->setCurrentText(currentCase.patient_gender);

    vtkNew<vtkDICOMImageReader> reader;
    QString dicomPath=currentCase.case_dir;
    reader->SetDirectoryName(dicomPath.toStdString().c_str());
    reader->Update();
}

void caseinfodialog::on_pushButton_close_clicked()
{
    close();
}
