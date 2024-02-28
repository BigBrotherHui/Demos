#ifndef MAINTOOLBAR_H
#define MAINTOOLBAR_H

#include "ufunbase.h"

namespace Ui {
class MainToolBar;
}

class MainToolBar : public uFunBase
{
    Q_OBJECT
public:
    static int typeId;

public:
    Q_INVOKABLE explicit MainToolBar(QWidget *parent = nullptr);
    ~MainToolBar();

private slots:

private:
    Ui::MainToolBar *ui;
};

#endif // MAINTOOLBAR_H
