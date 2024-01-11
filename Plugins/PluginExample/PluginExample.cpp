#include "PluginExample.h"
#include "widget.h"
#include <iostream>
#include "mesh_processing.h"
#include "ProjectionWidget.h"
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
	/*Eigen::MatrixX3d points;
	points.resize(8, 3);
	points.row(0) << 1, 0, 0;
	points.row(1) << 11, 0, 0;
	points.row(2) << 1, 10, 0;
	points.row(3) << 11, 10, 0;
	points.row(4) << 1, 0, 10;
	points.row(5) << 11, 0, 10;
	points.row(6) << 1, 10, 10;
	points.row(7) << 11, 10, 10;
	static_cast<ProjectionWidget*>(m_widget)->drawRect(points);*/
	return m_widget;
}
