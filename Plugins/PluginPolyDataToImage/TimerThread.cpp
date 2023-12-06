#include "TimerThread.h"
#include <QDebug>
TimerThread::TimerThread(QObject* parent) : QThread(parent)
{
}

TimerThread::~TimerThread()
{
	m_timer->stop();
	m_timer->deleteLater();
}

void TimerThread::run()
{
	m_timer = new QTimer;
	m_timer->setInterval(100);
	m_timer->start();
}

void TimerThread::onCreateTimer()
{
	m_timer = new QTimer();
	m_timer->setInterval(100);
	connect(m_timer, SIGNAL(timeout()), this, SLOT(onTimeout()));
	m_timer->start();
}

void TimerThread::onTimeout()
{
	qDebug() << "timerthread:"<<QThread::currentThreadId();
	emit signal_timeout();
}