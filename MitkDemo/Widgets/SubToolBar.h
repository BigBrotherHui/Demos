
#ifndef SubToolBar_H
#define SubToolBar_H

#include <QWidget>
#include "ufunbase.h"
#include <QtCore/qobjectdefs.h>

/**
  \brief SubToolBar

   IP-002 THA 主菜单

  \sa uFunBase
  \ingroup SubToolBar
*/

class QWidget;

namespace Ui
{
class SubToolBar;
}

class SubToolBar : public uFunBase
{
    Q_OBJECT

public:
    static int typeId;
    Q_INVOKABLE explicit SubToolBar(QWidget *parent = nullptr);
    ~SubToolBar();
    void f_Refresh() override;
    void f_Init() override;
private slots:
    void on_pushButton_open_clicked();

private:
    Ui::SubToolBar *ui;
};

#endif // SubToolBar_H
