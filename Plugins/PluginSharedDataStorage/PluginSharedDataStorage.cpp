#include "mainwindow.h"
#include "PluginSharedDataStorage.h"

int PluginSharedDataStorage::typeId = qRegisterMetaType<PluginSharedDataStorage*>("PluginSharedDataStorage");
PluginSharedDataStorage::PluginSharedDataStorage()
{
}

PluginSharedDataStorage::~PluginSharedDataStorage()
{
}

QWidget* PluginSharedDataStorage::createWidget()
{
	if(!m_widget)
		m_widget = new MainWindow;
	return m_widget;
}
