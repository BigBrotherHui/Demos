#include "PluginExample.h"
#include "widget.h"
#include <iostream>
Plugin::Plugin()
{
}

Plugin::~Plugin()
{
	if(m_widget)
	{
		delete m_widget;
		m_widget = nullptr;
	}
}

QWidget* Plugin::createWidget()
{
	if(!m_widget)
		m_widget = new Widget;
	return m_widget;
}
