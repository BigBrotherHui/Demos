
#ifndef StdMultiWidget_H
#define StdMultiWidget_H

#include <QWidget>
#include <QVBoxLayout>
#include "ufunbase.h"

/**
  \brief THARoboticsRegistrationsWidget

   图像显示MPR视图

  \sa uFunBase
  \ingroup THARoboticsRegistrationsWidget
*/

namespace Ui
{
class StdMultiWidget;
}
class QmitkLevelWindowWidget;
class QmitkStdMultiWidget;
class StdMultiWidget : public uFunBase
{
    Q_OBJECT

public:
    static int typeId;

    Q_INVOKABLE explicit StdMultiWidget(QWidget *parent = nullptr);
    ~StdMultiWidget();

    void f_Init() ;

    void f_Destroy();

    void f_Close();

    void f_Reset();

    int f_BeforeClose() override;

    QmitkStdMultiWidget *GetMultiWidget();
private:
    Ui::StdMultiWidget *ui;

    QVBoxLayout *vlayout;

    QWidget *viewParent;

    QHBoxLayout *hlayout ;

    std::unique_ptr<QmitkStdMultiWidget> multiWidget{nullptr};

    std::unique_ptr<QmitkLevelWindowWidget> levelWindowWidget;
};

#endif // StdMultiWidget_H
