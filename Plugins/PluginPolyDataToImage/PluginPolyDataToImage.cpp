#include <QDir>
#include <QDebug>
#include <QPluginLoader>
#include "PluginInterface.h"
#include "PluginPolyDataToImage.h"
#include "ReamWidget.h"
PluginPolyDataToImage::PluginPolyDataToImage()
{
}

PluginPolyDataToImage::~PluginPolyDataToImage()
{
	if (m_widget)
	{
		delete m_widget;
		m_widget = nullptr;
	}
}

QWidget* PluginPolyDataToImage::createWidget()
{
    if (!m_widget)
        m_widget = new ReamWidget;
    return m_widget;
}