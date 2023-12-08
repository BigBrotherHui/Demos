#include "mainwindow.h"
#include "PluginSharedDataStorage.h"

PluginSharedDataStorage::PluginSharedDataStorage()
{
}

PluginSharedDataStorage::~PluginSharedDataStorage()
{
	if (m_widget)
	{
		delete m_widget;
		m_widget = nullptr;
	}
}

QWidget* PluginSharedDataStorage::createWidget()
{
	if(!m_widget)
		m_widget = new MainWindow;
	return m_widget;
}
