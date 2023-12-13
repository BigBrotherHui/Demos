#include "PluginExample.h"
#include "widget.h"
#include <iostream>
int PluginExample::typeId = qRegisterMetaType<PluginExample*>("PluginExample");
PluginExample::PluginExample()
{
}

PluginExample::~PluginExample()
{
}

QWidget* PluginExample::createWidget()
{
	if(!m_widget)
		m_widget = new Widget;
	return m_widget;
}
