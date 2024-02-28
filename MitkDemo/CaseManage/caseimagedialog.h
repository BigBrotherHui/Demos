#ifndef CASEIMAGEDIALOG_H
#define CASEIMAGEDIALOG_H

#include <QDialog>
#include "vtkSmartPointer.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkImageViewer2.h"
#include <QPointer>
#include "seriesselectdialog.h"
namespace Ui {
class CaseImageDialog;
}

class CaseImageDialog : public QDialog
{
    Q_OBJECT

public:
    explicit CaseImageDialog(QWidget *parent = nullptr);
    ~CaseImageDialog();
    void SetCurrentCaseId(int pCaseId) { mCaseId=pCaseId; }
private slots:
    void on_pushButton_image_clicked();

    void on_pushButton_close_clicked();

    void on_verticalSlider_slice_valueChanged(int pValue);

    void on_pushButton_import_clicked();


private:
    Ui::CaseImageDialog *ui;
    QPointer<SeriesSelectDialog> dialog{nullptr};
    vtkSmartPointer< vtkImageViewer2> mImageView ;
    int mCaseId{-1};
};

#endif // CASEIMAGEDIALOG_H
