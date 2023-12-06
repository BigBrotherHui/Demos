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
}

QWidget* PluginPolyDataToImage::createWidget()
{
    if (!m_widget)
        m_widget = new ReamWidget;
    return m_widget;
}