
#pragma once

#include <QWidget>
#include "ufunbase.h"
#include "ufunction.h"
#include "ustatus.h"
#include <stack>
#include <vtkImageData.h>
#include <vtkDICOMImageReader.h>
class vtkImageViewer2;
class QPushButton;
class QTableWidget;
namespace Ui {
class CaseManageControl;
}
class CaseManageControl : public uFunBase
{
    Q_OBJECT

public:
    static int typeId;
    struct RenderInfo
    {
    public:
        QString imageDir;
        int slice;
        double scale{-1};
        static vtkSmartPointer<vtkDICOMImageReader> reader;
        RenderInfo(const RenderInfo &info)
        {
            this->slice=info.slice;
            this->scale=info.scale;
            this->imageDir=info.imageDir;
        }
        RenderInfo(){}
    };
    Q_INVOKABLE explicit CaseManageControl(QWidget *parent = nullptr);
    ~CaseManageControl();
protected:
    void ShowDicom();
    Q_INVOKABLE void enableButtons(bool pIsEnable);
private slots:
    void slot_updateDicom(QTableWidget* w, QTableWidgetItem* cur, QTableWidgetItem* prev);
    void on_cb_new_clicked();

    void on_cb_modi_clicked();

    void on_cb_open_clicked();

    void on_pushButton_search_clicked();

    void on_lineEdit_search_returnPressed();

    void on_cb_del_clicked();


    void on_verticalSlider_slice_valueChanged(int value);

    void on_checkBox_previewDicom_clicked();

    void on_cb_close_clicked();

    void on_pushButton_image_clicked();

private:
    Ui::CaseManageControl*ui;
    std::stack<RenderInfo> mRenderInfoVec;
    vtkSmartPointer<vtkImageViewer2> mImageViewer{nullptr};
    int mLastSelectedRow{-1};
    QTimer *mRenderTimer{nullptr};
    QString mCurrentImageDir;
    QList<QPushButton *> mEnableButtons;
};