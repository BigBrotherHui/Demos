#include "PluginExample.h"
#include "widget.h"
Plugin::Plugin()
{
}

QWidget* Plugin::createWidget()
{
	if(!m_widget)
		m_widget = new Widget;
	return m_widget;
}
