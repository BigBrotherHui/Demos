#ifndef TIMERTHREAD_H
#define TIMERTHREAD_H
#include <QThread>
#include <QTimer>

class TimerThread : public QThread
{
    Q_OBJECT
public:
    TimerThread(QObject *parent=nullptr);
    ~TimerThread();
    virtual void run()override;
public slots:
    void onCreateTimer();
    void onTimeout();
signals:
    void signal_timeout();
private:
    QTimer* m_timer;
};

#endif // MYTHREAD_H
