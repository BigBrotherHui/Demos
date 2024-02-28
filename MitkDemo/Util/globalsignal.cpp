#include "globalsignal.h"
#include <QTableWidget>i
GlobalSignal::GlobalSignal(QObject *parent) : QObject(parent)
{

}

void GlobalSignal::slot_updateDicom(QTableWidgetItem *cur, QTableWidgetItem *prev)
{
    QTableWidget *w=dynamic_cast<QTableWidget *>(sender());
    if(!w)
        return;
    emit signal_updateDicom(w,cur,prev);
}
