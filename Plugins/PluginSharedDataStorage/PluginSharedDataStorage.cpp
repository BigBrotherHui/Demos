#include "mainwindow.h"
#include "PluginSharedDataStorage.h"

PluginSharedDataStorage::PluginSharedDataStorage()
{
}

QWidget* PluginSharedDataStorage::createWidget()
{
	if(!m_widget)
		m_widget = new MainWindow;
	return m_widget;
}
