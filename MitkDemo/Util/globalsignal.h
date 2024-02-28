#ifndef GLOBALSIGNAL_H
#define GLOBALSIGNAL_H

#include <QObject>
class QTableWidgetItem;
class QTableWidget;
class GlobalSignal : public QObject
{
    Q_OBJECT
public:
    static GlobalSignal& instance()
    {
        static GlobalSignal signal;
        return signal;
    }
    void slot_updateDicom(QTableWidgetItem *cur,QTableWidgetItem *prev);
signals:
    void signal_updateDicom(QTableWidget *w,QTableWidgetItem *cur,QTableWidgetItem *prev);
protected:
    GlobalSignal(QObject *parent=nullptr);
};

#endif // GLOBALSIGNAL_H
