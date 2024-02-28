#ifndef UBLANK_H
#define UBLANK_H

#include <QWidget>
#include "ufunbase.h"

namespace Ui {
class uBlank;
}

class uBlank : public uFunBase
{
    Q_OBJECT

public:
    static int typeId;

    Q_INVOKABLE explicit uBlank(QWidget *parent = nullptr);
    ~uBlank();

private:
    Ui::uBlank *ui;
};

#endif // UBLANK_H
