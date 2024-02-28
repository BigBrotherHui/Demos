#ifndef SERIESSELECTDIALOG_H
#define SERIESSELECTDIALOG_H

#include <QDialog>
#include <vtkDICOMDirectory.h>
#include <vtkDICOMReader.h>
#include <vtkNew.h>
#include "seriesinfo.h"

class CaseImageDialog;
class DicomImageUtil;
namespace Ui {
class SeriesSelectDialog;
}

class SeriesSelectDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SeriesSelectDialog(QWidget *parent = nullptr);
    ~SeriesSelectDialog();
    vtkStringArray *GetFileNamesForSeries();
Q_SIGNALS:
    /**
     * @brief 加载图像
     * @param pFiles 要加载的图像列表
     * @param pSeriesInfo 要加载的图像序列信息
     */
    void SignalLoadImage(vtkStringArray* pFiles, SeriesTuple pSeriesInfo);
public Q_SLOTS:
    /**
     * @brief 加载序列列表
     * @param folder 图像所在文件夹
     */
    void loadSeries(QString folder);

private slots:
    void on_pushButton_cancel_clicked();

    void on_pushButton_ok_clicked();

    void on_tableWidget_series_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn);

    void on_tableWidget_series_doubleClicked(const QModelIndex &index);

    void on_pushButton_close_clicked();

private:
    Ui::SeriesSelectDialog *ui;
    CaseImageDialog* mCaseImageDialog;
    DicomImageUtil* mDicomImageUtil;
    int COLUMNS_COUNT = 4;
    int COLUMN_PREVIEW = 0;
    int COLUMN_IMAGE_COUNT = 1;
    int COLUMN_INFO = 2;
    int COLUMN_CHOOSE = 3;
    int mSelectedSeries = -1;
    QList<SeriesTuple>* mSeriesList;
    vtkNew<vtkDICOMDirectory> mDicomDirectory;
    vtkNew<vtkDICOMReader> mDICOMReader ;
    vtkNew<vtkStringArray> mStringArray;
    /**
     * @brief 获取预览图标
     * @param image 原始图像
     * @return 预览图标
     */
    QIcon getQIcon(QImage image);
    /**
     * @brief 获取预览图像
     * @param imageData 原始dicom图像数据
     * @return 预览图像
     */
    QImage getQImage(vtkImageData *imageData);
    /**
     * @brief 转换dicom图像
     * @param inPtr 输入数据
     * @param outPtr 输出数据
     * @param ncols 列数
     * @param nrows 行数
     * @param numComp通道数
     */
    void ConvertVtkImageToRGBA(const unsigned char *inPtr, unsigned char *outPtr, int ncols, int nrows, int numComp);
    /**
     * @brief 设置行选择状态
     * @param rowIndex 行号
     * @param selected 是否选中
     */
    void setRowSelect(int rowIndex, bool selected);
};

#endif // SERIESSELECTDIALOG_H
