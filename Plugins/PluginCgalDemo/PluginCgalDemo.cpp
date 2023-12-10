#include <QDir>
#include <QDebug>
#include <QPluginLoader>
#include "PluginInterface.h"
#include "PluginCgalDemo.h"
#include "Widget.h"
PluginCgalDemo::PluginCgalDemo()
{
}

PluginCgalDemo::~PluginCgalDemo()
{
	if (m_widget)
	{
		delete m_widget;
		m_widget = nullptr;
	}
}

QWidget* PluginCgalDemo::createWidget()
{
    if (!m_widget)
        m_widget = new Widget(nullptr);
    return m_widget;
}