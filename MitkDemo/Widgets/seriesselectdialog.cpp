#pragma execution_character_set("utf-8")

#include "seriesselectdialog.h"
#include "ui_seriesselectdialog.h"

#include "QMessageBox.h"
#include <QJsonObject>
#include <QTableWidgetItem>
#include <QDebug>

#include "vtkStringArray.h"
#include "vtkImageData.h"
#include "vtkImageShiftScale.h"
#include "progressutil.h"
#include "dicomimageutil.h"
#include "udicomtag.h"
#include "caseimagedialog.h"
#include <QStyledItemDelegate>
#include <QPainter>
class Delegate2 : public QStyledItemDelegate
{
public:
    Delegate2(QObject *parent = nullptr) : QStyledItemDelegate(parent){}
    void paint(QPainter* painter, const QStyleOptionViewItem& option,
          const QModelIndex& index) const
    {
        QStyleOptionViewItem opt(option);
        if (opt.state & QStyle::State_HasFocus)
            opt.state = opt.state & ~QStyle::State_HasFocus;
        if(index.column()==3)
        {
            QIcon icon(index.data(Qt::DecorationRole).value<QIcon>());
            QPixmap pixmap=icon.pixmap(icon.availableSizes().first());
            if(pixmap.width()>19)
                pixmap=pixmap.scaled(17,17);
            QBrush brush(opt.palette.background());
            if(opt.state & QStyle::State_Selected)
                brush=QBrush("#1167CF");
            painter->setBrush(brush);
            painter->setPen(Qt::transparent);
            painter->drawRect(opt.rect);
            qApp->style()->drawItemPixmap(painter,opt.rect,Qt::AlignCenter,pixmap);
            QPen p(QColor("#303439"),1);
            painter->setPen(p);
            painter->drawLine(opt.rect.bottomLeft(),opt.rect.bottomRight());
            return;
        }
        return QStyledItemDelegate::paint(painter, opt, index);
    }
};

SeriesSelectDialog::SeriesSelectDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SeriesSelectDialog)
{
    ui->setupUi(this);
    mCaseImageDialog = static_cast<CaseImageDialog*>(parent);
    mDicomImageUtil = new DicomImageUtil;
    ui->tableWidget_series->horizontalHeader()->setDefaultAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    ui->tableWidget_series->setEditTriggers(QAbstractItemView::NoEditTriggers);
    ui->tableWidget_series->horizontalHeader()->setStretchLastSection(true);
    ui->tableWidget_series->setSelectionBehavior(QAbstractItemView::SelectRows);
    ui->tableWidget_series->horizontalHeader()->setHighlightSections(false);
    ui->tableWidget_series->horizontalHeader()->setFixedHeight(40);
    ui->tableWidget_series->verticalHeader()->setDefaultSectionSize(80);
    QStringList columnLabels = {tr("预览"), tr("图像数量"), tr("信息"), tr("选择")};
    ui->tableWidget_series->setColumnCount(columnLabels.size());
    ui->tableWidget_series->setHorizontalHeaderLabels(columnLabels);
    ui->tableWidget_series->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
    ui->tableWidget_series->setColumnWidth(COLUMN_PREVIEW, 80);
    ui->tableWidget_series->setColumnWidth(COLUMN_IMAGE_COUNT, 100);
    ui->tableWidget_series->setColumnWidth(COLUMN_INFO, 200);
    ui->tableWidget_series->setColumnWidth(COLUMN_CHOOSE, 80);
    ui->tableWidget_series->setIconSize(QSize(64,64));
    ui->tableWidget_series->verticalHeader()->setVisible(true);
    ui->tableWidget_series->verticalHeader()->setFixedWidth(44);
    ui->tableWidget_series->verticalHeader()->setDefaultAlignment(Qt::AlignCenter);
    ui->tableWidget_series->setItemDelegate(new Delegate2(this));

    ui->pushButton_ok->setFixedSize(96,34);
    ui->pushButton_cancel->setFixedSize(96,34);
}

SeriesSelectDialog::~SeriesSelectDialog()
{
    disconnect(this, 0, 0, 0);
    delete ui;
}

vtkStringArray *SeriesSelectDialog::GetFileNamesForSeries()
{
    return mStringArray;
}


void SeriesSelectDialog::on_pushButton_cancel_clicked()
{
    this->reject();
}

void SeriesSelectDialog::on_pushButton_ok_clicked()
{
    if(mSelectedSeries<0){
        QMessageBox::information(this, tr("提示"), tr("请选择序列"));
        return;
    }
    SeriesTuple seriesInfo = mSeriesList->at(mSelectedSeries);
    mStringArray->DeepCopy(mDicomDirectory->GetFileNamesForSeries(mSelectedSeries));
    if(!mDicomImageUtil->CheckSeries(mStringArray))
    {
        QMessageBox::information(this, tr("提示"), tr("选择的序列有丢帧，请联系图像提供方"));
        return;
    }
    emit SignalLoadImage(mStringArray, seriesInfo);
    this->accept();
}

void SeriesSelectDialog::on_tableWidget_series_currentCellChanged(int currentRow, int , int previousRow, int )
{
    if(currentRow==previousRow){
        return;
    }
    setRowSelect(previousRow, false);
    setRowSelect(currentRow, true);
    mSelectedSeries = currentRow;
}

void SeriesSelectDialog::loadSeries(QString folder){
    qDebug()<<"SeriesSelectDialog: load series";
    if(folder.isEmpty()){
        return;
    }
    int step = 0, stepScan = 40;
    ProgressUtil progressUtil;
    progressUtil.Start(tr("检索图像序列中..."));
    ui->tableWidget_series->setRowCount(0);
    mSeriesList = new QList<SeriesTuple>();
    mDicomDirectory->SetDirectoryName(folder.toStdString().c_str());
    //设置自动递归查询文件下的dicom数据  3层
    mDicomDirectory->SetScanDepth(3);
    mDicomDirectory->Update();
    step += stepScan;
    progressUtil.Step(step);
    //序列数
    int seriesCount = mDicomDirectory->GetNumberOfSeries();
    if(seriesCount<=0)
    {
        progressUtil.Finish();
        return;
    }
    int stepLast = 100 - step;
    for(int seriesIndex = 0,rowIndex = 0;seriesIndex<seriesCount; seriesIndex++,rowIndex++)
    {
        SeriesTuple seriesInfo = {};
        ui->tableWidget_series->insertRow(rowIndex);
        vtkStringArray* stringArray = mDicomDirectory->GetFileNamesForSeries(seriesIndex);
        //切片数量
        seriesInfo.ImageCount = stringArray->GetNumberOfValues();
        uDicomTag* dicomTag = new uDicomTag();
        dicomTag->f_OpenFileArray(stringArray);
        seriesInfo.AccessNo = QString::fromStdString(dicomTag->f_GetAccessionNumber());
        seriesInfo.SeriesNo = QString::fromStdString(dicomTag->f_GetSeriesNumber());
        seriesInfo.WindowWidth = atoi(dicomTag->f_GetWindowWidth().c_str());
        seriesInfo.WindowCenter = atoi(dicomTag->f_GetWindowCenter().c_str());
        mSeriesList->append(seriesInfo);
        for (int columnIndex = 0; columnIndex < COLUMNS_COUNT; ++columnIndex)
        {
            QTableWidgetItem* item = new QTableWidgetItem();
            if(columnIndex==COLUMN_IMAGE_COUNT){
                item->setText(QString::number(seriesInfo.ImageCount));
            }
            else if(columnIndex==COLUMN_INFO){
                //获取信息
                QString patientName = QString::fromStdString(dicomTag->f_GetPatientName());
                QString accessNo = QString::fromStdString(dicomTag->f_GetAccessionNumber());
                QString seriesNo = QString::fromStdString(dicomTag->f_GetSeriesNumber());
                QString info;
                info.append(tr("PatientName:")).append(patientName).append(tr("\r\nImageNo:")).append(accessNo).append(tr("\r\nSeriesNo:")).append(seriesNo);
                item->setText(info);
            }
            else if(columnIndex==COLUMN_PREVIEW){
                item->setIcon(getQIcon(getQImage(dicomTag->f_GetImageData())));
            }
            item->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
            ui->tableWidget_series->setItem(rowIndex, columnIndex, item);
        }
        setRowSelect(rowIndex, false);
        step += stepLast/seriesCount;
    }
    progressUtil.Finish();
}

void SeriesSelectDialog::ConvertVtkImageToRGBA(const unsigned char *inPtr, unsigned char *outPtr, int ncols, int nrows, int numComp)
{
    unsigned char alpha = 255;

    int inIncX = numComp;

    // 遍历数据，并复制纹理
    if (numComp == 1)
    {
        for (int idy = 0; idy < nrows; idy++)
        {
            for (int idx = 0; idx < ncols; idx++)
            {
                // 引导操作
                int val = inPtr[0];

                unsigned char cval = static_cast<unsigned char>(val + 0.5);
                outPtr[0] = cval;
                outPtr[1] = cval;
                outPtr[2] = cval;
                outPtr[3] = alpha;
                outPtr += 4;
                inPtr += inIncX;
            }
        }
    }
    else if (numComp == 2)
    {
        for (int idy = 0; idy < nrows; idy++)
        {
            for (int idx = 0; idx < ncols; idx++)
            {
                // 引导操作
                int val = inPtr[0];
                unsigned char cval = static_cast<unsigned char>(val + 0.5);
                val = inPtr[1];
                unsigned char aval = static_cast<unsigned char>(val + 0.5);
                outPtr[0] = cval;
                outPtr[1] = cval;
                outPtr[2] = cval;
                outPtr[3] = aval;
                outPtr += 4;
                inPtr += inIncX;
            }
        }
    }
    else if (numComp == 3)
    {
        for (int idy = 0; idy < nrows; idy++)
        {
            for (int idx = 0; idx < ncols; idx++)
            {
                // 引导操作
                int r = inPtr[0];
                int g = inPtr[1];
                int b = inPtr[2];
                outPtr[0] = static_cast<unsigned char>(r + 0.5);
                outPtr[1] = static_cast<unsigned char>(g + 0.5);
                outPtr[2] = static_cast<unsigned char>(b + 0.5);
                outPtr[3] = alpha;
                outPtr += 4;
                inPtr += inIncX;
            }
        }
    }
    else // if (numComp == 4)
    {
        for (int idy = 0; idy < nrows; idy++)
        {
            for (int idx = 0; idx < ncols; idx++)
            {
                // Pixel operation
                int r = inPtr[0];
                int g = inPtr[1];
                int b = inPtr[2];
                int a = inPtr[3];
                outPtr[0] = static_cast<unsigned char>(r + 0.5);
                outPtr[1] = static_cast<unsigned char>(g + 0.5);
                outPtr[2] = static_cast<unsigned char>(b + 0.5);
                outPtr[3] = static_cast<unsigned char>(a + 0.5);
                outPtr += 4;
                inPtr += inIncX;
            }
        }
    }
}

QImage SeriesSelectDialog::getQImage(vtkImageData *imageData)
{
     // 图像列数（宽）
    int ncols = imageData->GetDimensions()[0];
    // 图像行数（高）
    int nrows = imageData->GetDimensions()[1];
    // 获取原图像数据的灰度值范围
    // 如：CT值的单位是Hounsfield，简称为Hu，范围是-1024-3071。用于衡量人体组织对X射线的吸收率，设定水的吸收率为0Hu。
    double range[2];
    imageData->GetScalarRange(range);
    // 将原图像灰度值调节到0~255范围内，输出为unsigned char（0~255）类型数据（即调节窗宽窗位）
    vtkNew<vtkImageShiftScale> shifter ;
    //减去最小值
    shifter->SetShift(-1.0*range[0]);
    //调节比例，灰度换算公式：255*X/(range[1]-range[0])
    shifter->SetScale(255.0 / (range[1] - range[0]));
    shifter->SetOutputScalarTypeToUnsignedChar();
    shifter->SetInputData(imageData);
    shifter->Update();
    // 获取原图像component size值，即通道数
    int numComp = shifter->GetOutput()->GetScalarSize();
    // 申请内存
    // 输入，从图像中获取数据指针
    unsigned char *inPtr = static_cast<unsigned char*>(shifter->GetOutput()->GetScalarPointer());
    // 输出，RGBA共4通道，所以这里乘4
    unsigned char *outPtr = new unsigned char[nrows * ncols * 4];
    ConvertVtkImageToRGBA(inPtr, outPtr, ncols, nrows, numComp);
    // 创建QImage，理论上可以用outPtr创建任意其它类型图像
    QImage image = QImage(outPtr, ncols, nrows, QImage::Format_RGB32);
    return image;
}

QIcon SeriesSelectDialog::getQIcon(QImage image)
{
    return QIcon(QPixmap::fromImage(image).scaled(QSize(64,64),Qt::IgnoreAspectRatio,Qt::SmoothTransformation));
}

void SeriesSelectDialog::setRowSelect(int rowIndex, bool selected){
    if (rowIndex < 0) {
        return;
    }
    if (rowIndex >= this->ui->tableWidget_series->rowCount()) {
        return;
    }
    QTableWidgetItem* ckBoxItem =
        static_cast<QTableWidgetItem *>(this->ui->tableWidget_series->item(rowIndex, COLUMN_CHOOSE));
    if (ckBoxItem == nullptr) {
        ckBoxItem = new QTableWidgetItem();
        this->ui->tableWidget_series->setItem(rowIndex, COLUMN_CHOOSE, ckBoxItem);
    }
    ckBoxItem->setCheckState(selected? Qt::Checked:Qt::Unchecked);
}

void SeriesSelectDialog::on_tableWidget_series_doubleClicked(const QModelIndex &)
{
    on_pushButton_ok_clicked();
}

void SeriesSelectDialog::on_pushButton_close_clicked()
{
    close();
}
