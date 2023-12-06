#include "PluginExample.h"
#include "widget.h"
Plugin::Plugin()
{
}

QWidget* Plugin::widget()
{
	if(!m_widget)
		m_widget = new Widget;
	return m_widget;
}
