#include <QDir>
#include <QDebug>
#include <QPluginLoader>
#include "PluginInterface.h"
#include "PluginPolyDataToImage.h"
#include "ReamWidgetTest.h"
int PluginPolyDataToImage::typeId = qRegisterMetaType<PluginPolyDataToImage*>("PluginPolyDataToImage");
PluginPolyDataToImage::PluginPolyDataToImage()
{
}

PluginPolyDataToImage::~PluginPolyDataToImage()
{
}

QWidget* PluginPolyDataToImage::createWidget()
{
    if (!m_widget)
        m_widget = new ReamWidgetTest;
    return m_widget;
}