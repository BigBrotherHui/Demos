#include "PluginExample.h"
#include "widget.h"
#include <iostream>
#include "mesh_processing.h"
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
		m_widget = new MeshProcessing;
	return m_widget;
}
